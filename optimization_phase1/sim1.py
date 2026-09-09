import os
import math

# =================================================================
# --- 1. HELPER FUNCTIONS (SI UNITS: m, N, Pa, kg)
# =================================================================
def generate_dd_stack(n_reps, phi, psi):
    """Generates a Double-Double laminate stack sequence: [phi, -psi, -phi, psi] x n_reps."""
    if n_reps <= 0:
        return []
    block = [phi, -psi, -phi, psi]
    return block * int(n_reps)

def calculate_D_matrix(E1, E2, nu12, G12, stack, ply_thickness=0.000127):
    """
    Calculates the bending stiffness matrix (D-matrix) for a laminate stack in N*m.
    
    Parameters:
        E1, E2        : Longitudinal and transverse Young's moduli (Pa)
        nu12          : Major Poisson's ratio
        G12           : In-plane shear modulus (Pa)
        stack         : List of ply orientation angles in degrees
        ply_thickness : Nominal thickness of a single ply (m)
    """
    # Calculate minor Poisson's ratio and reduced stiffness matrix (Q-matrix)
    nu21 = (E2 / E1) * nu12
    Q11 = E1 / (1.0 - nu12 * nu21)
    Q22 = E2 / (1.0 - nu12 * nu21)
    Q12 = (nu12 * E2) / (1.0 - nu12 * nu21)
    Q66 = G12
    
    D11 = D22 = D12 = D66 = 0.0
    total_thickness = len(stack) * ply_thickness
    z_current = -total_thickness / 2.0
    
    # Integrate transformed Q-matrix (Q-bar) across ply coordinates along z-axis
    for angle in stack:
        theta = math.radians(angle)
        m, n = math.cos(theta), math.sin(theta)
        m2, n2, m4, n4, m2n2 = m**2, n**2, m**4, n**4, (m**2)*(n**2)
        
        # Transformed reduced stiffness components
        Q11_bar = Q11*m4 + 2.0*(Q12 + 2.0*Q66)*m2n2 + Q22*n4
        Q22_bar = Q11*n4 + 2.0*(Q12 + 2.0*Q66)*m2n2 + Q22*m4
        Q12_bar = (Q11 + Q22 - 4.0*Q66)*m2n2 + Q12*(m4 + n4)
        Q66_bar = (Q11 + Q22 - 2.0*Q12 - 2.0*Q66)*m2n2 + Q66*(m4 + n4)
        
        # Ply lower and upper z-coordinates
        z_k_minus_1 = z_current
        z_k = z_current + ply_thickness
        z_factor = (z_k**3 - z_k_minus_1**3) / 3.0
        
        # Accumulate bending stiffness contributions
        D11 += Q11_bar * z_factor
        D22 += Q22_bar * z_factor
        D12 += Q12_bar * z_factor
        D66 += Q66_bar * z_factor
        z_current = z_k
        
    return D11, D22, D12, D66

def calculate_analytical_buckling(D11, D22, D12, D66, a, b, Nx, Ny, Nxy):
    """
    Calculates the critical buckling load factor (Lambda_c) for a rectangular panel 
    under combined biaxial and shear in-plane loads.
    
    Parameters:
        D11, D22, D12, D66 : Bending stiffness components (N*m)
        a, b               : Panel length and width dimensions (m)
        Nx, Ny, Nxy        : Applied in-plane resultant forces per unit length (N/m)
    """
    lambda_n_min = float('inf')
    
    # 1. Uniaxial / Biaxial Buckling Mode Iteration over half-wave numbers (m, n)
    for m in range(1, 15):
        for n in range(1, 5):
            num = (math.pi**2) * (D11 * (m/a)**4 + 2.0*(D12 + 2.0*D66) * ((m/a)**2) * ((n/b)**2) + D22 * (n/b)**4)
            den = ((m/a)**2) * Nx + ((n/b)**2) * Ny
            if den > 0:
                lambda_n = num / den
                if lambda_n < lambda_n_min:
                    lambda_n_min = lambda_n

    if lambda_n_min == float('inf'):
        lambda_n_min = 1e6 

    # 2. Pure Shear Buckling Factor (Lambda_s) Estimation
    Nxy_abs = abs(Nxy)
    if Nxy_abs < 1e-6:
        lambda_s = 1e6
    else:
        gamma = math.sqrt(D11 * D22) / (D12 + 2.0*D66)
        beta_1 = 10.0
        if gamma >= 1.0:
            lambda_s = (4.0 * beta_1 * (D11 * D22**3)**0.25) / (b**2 * Nxy_abs)
        else:
            lambda_s = (4.0 * beta_1 * math.sqrt(D22 * (D12 + 2.0*D66))) / (b**2 * Nxy_abs)

    # 3. Combined Load Interaction Equation (Biaxial + Shear)
    A = 1.0 / (lambda_s**2)
    B = 1.0 / lambda_n_min
    if A > 1e-12:
        lambda_c = (-B + math.sqrt(B**2 + 4.0*A)) / (2.0*A)
    else:
        lambda_c = lambda_n_min
        
    return lambda_c

# =================================================================
# --- 2. HYPERMESH MODEL GENERATOR (MID-THICKNESS TIP LOADS)
# =================================================================
def generate_wing_box_model(
    dd_angles=[35, 60], 
    lower_reps={0: 6, 1: 5, 2: 4}, 
    upper_reps={0: 6, 1: 5, 2: 4}, 
    spars_material="aluminum", 
    spars_plies_override=None
):
    """
    Generates a Tcl script for Altair HyperMesh to set up materials, PCOMP/PSHELL properties,
    mesh classification, boundary conditions, tip loads, and static/buckling subcases.
    """
    OUTPUT_FILE = "optimum.tcl"
    DEFAULT_MESH_SIZE = 0.0635   # 63.5 mm mesh resolution
    MESH_TYPE = 2
    PLY_THICKNESS_SI = 0.000127  # 0.127 mm per ply
    
    # GEOMETRY PARAMETERS (SI UNITS)
    Z_LOWER = 0.0                  # 0 m
    Z_UPPER = 0.3810               # 15 in (0.3810 m)
    Z_MID   = Z_UPPER / 2.0        # 7.5 in (0.1905 m)
    
    # Material Definitions
    materials = {
        "Composite": { 
            "id": 1, "type": "MAT8", "name": "T300_5208", 
            "E1": 127.55e9, "E2": 13.03e9, "NU12": 0.3, "G12": 6.412e9, "RHO": 1580.0 
        },
        "Aluminum":  { 
            "id": 2, "type": "MAT1", "name": "Al_2024_T3", 
            "E": 68.95e9, "NU": 0.3, "RHO": 2770.0 
        }
    }

    zones = []
    # Grid indexing map for skin panel layout (3x3 grid for upper & lower covers)
    panel_grid_map = {
        1: (0, 0), 2: (0, 1), 3: (0, 2),  4: (1, 0), 5: (1, 1), 6: (1, 2),  7: (2, 0), 8: (2, 1), 9: (2, 2), 
        10: (0, 0), 11: (0, 1), 12: (0, 2), 13: (1, 0), 14: (1, 1), 15: (1, 2), 16: (2, 0), 17: (2, 1), 18: (2, 2)
    }

    # Lower Panels (IDs 1-9)
    for panel_id in range(1, 10):
        col, row = panel_grid_map[panel_id] 
        reps = lower_reps[col] if isinstance(lower_reps, dict) else lower_reps
        final_stack = generate_dd_stack(reps, dd_angles[0], dd_angles[1])
        zones.append({
            "name": f"PCOMP_Low_Zone{col+1}", 
            "component": f"Lower_Panel_{panel_id}", 
            "stack": final_stack, 
            "thickness": PLY_THICKNESS_SI 
        })

    # Upper Panels (IDs 10-18)
    for panel_id in range(10, 19):
        col, row = panel_grid_map[panel_id]
        reps = upper_reps[col] if isinstance(upper_reps, dict) else upper_reps
        final_stack = generate_dd_stack(reps, dd_angles[0], dd_angles[1])
        zones.append({
            "name": f"PCOMP_Up_Zone{col+1}", 
            "component": f"Upper_Panel_{panel_id}", 
            "stack": final_stack, 
            "thickness": PLY_THICKNESS_SI 
        })

    print(f"Generating TCL script (Tip loads applied at Z_mid = {Z_MID} m)...")
    
    with open(OUTPUT_FILE, "w", encoding='utf-8') as f:
        f.write("# Generator Script v8 - Tip Loads applied at Mid-Thickness (Z_mid = 0.1905 m)\n")
        f.write("hm_answernext yes\n")
        f.write("hm_blockmessages 1\n\n")

        # 1. CLEANUP PREVIOUS ENTITIES
        f.write("catch { *clearmark loadsteps 1; *createmark loadsteps 1 all; *deletemark loadsteps 1 }\n")
        f.write("catch { *clearmark loadstepinputs 1; *createmark loadstepinputs 1 all; *deletemark loadstepinputs 1 }\n")
        f.write("catch { *clearmark loadcols 1; *createmark loadcols 1 all; *deletemark loadcols 1 }\n")
        f.write("catch { *clearmark props 1; *createmark props 1 all; *deletemark props 1 }\n")
        f.write("catch { *clearmark mats 1; *createmark mats 1 all; *deletemark mats 1 }\n\n")

        # 2. SURFACE MESHING AND EQUIVALENCING
        f.write("*createmark surfs 1 all\n")
        f.write("catch { *edgeequivalence surfs 1 0.00254 }\n")
        f.write(f"catch {{ *defaultremeshsurf 1 {DEFAULT_MESH_SIZE} {MESH_TYPE} 2 2 1 1 1 1 0 0 0 0 }}\n")
        f.write("*createmark nodes 1 all\n")
        f.write("if {[hm_marklength nodes 1] > 0} { catch { *equivalence nodes 1 0.0254 1 0 0 } }\n\n")

        # 3. INITIALIZE COMPONENTS
        for zone in zones: 
            f.write(f"if {{[hm_entityinfo exist comps \"{zone['component']}\" -byname] == 0}} {{ *createentity comps name=\"{zone['component']}\" }}\n")
            f.write(f"set list_{zone['component']} \"\"\n")
            
        f.write("if {[hm_entityinfo exist comps \"Spars_and_Ribs\" -byname] == 0} { *createentity comps name=\"Spars_and_Ribs\" }\n")
        f.write("set list_Spars_and_Ribs \"\"\n\n")

        # Element spatial classification: Z_LOWER = 0.0m, Z_UPPER = 0.3810m
        f.write("*createmark elems 1 all\n")
        f.write("set all_elems [hm_getmark elems 1]\n")
        f.write("foreach elem_id $all_elems {\n")
        f.write("    set n [hm_getvalue elems id=$elem_id dataname=nodes]\n")
        f.write("    set min_z 99999.0; set max_z -99999.0\n")
        f.write("    set sum_x 0.0; set sum_y 0.0; set sum_z 0.0\n")
        f.write("    set num_n [llength $n]\n")
        f.write("    foreach n_id $n {\n")
        f.write("        set nx [hm_getvalue nodes id=$n_id dataname=x]\n")
        f.write("        set ny [hm_getvalue nodes id=$n_id dataname=y]\n")
        f.write("        set nz [hm_getvalue nodes id=$n_id dataname=z]\n")
        f.write("        set sum_x [expr {$sum_x + $nx}]\n")
        f.write("        set sum_y [expr {$sum_y + $ny}]\n")
        f.write("        set sum_z [expr {$sum_z + $nz}]\n")
        f.write("        if {$nz < $min_z} { set min_z $nz }\n")
        f.write("        if {$nz > $max_z} { set max_z $nz }\n")
        f.write("    }\n")
        f.write("    set cx [expr {$sum_x / double($num_n)}]\n")
        f.write("    set cy [expr {$sum_y / double($num_n)}]\n")
        f.write("    set cz [expr {$sum_z / double($num_n)}]\n")
        f.write("    set delta_z [expr {$max_z - $min_z}]\n\n")
        
        # Planar filter distinguishing upper/lower skins from internal spars/ribs
        f.write("    set is_skin 0\n")
        f.write("    set is_upper 0\n")
        f.write("    if {$delta_z < 0.005} {\n")
        f.write(f"        if {{abs($cz - {Z_LOWER}) < 0.01}} {{\n")
        f.write("            set is_skin 1; set is_upper 0\n")
        f.write(f"        }} elseif {{abs($cz - {Z_UPPER}) < 0.01}} {{\n")
        f.write("            set is_skin 1; set is_upper 1\n")
        f.write("        }\n")
        f.write("    }\n\n")
        
        f.write("    if {$is_skin == 1} {\n")
        f.write("        set col -1; set row -1\n")
        f.write("        if {$cx >= -0.01 && $cx <= 1.1811} { set col 0 } elseif {$cx > 1.1811 && $cx <= 2.3622} { set col 1 } elseif {$cx > 2.3622 && $cx <= 3.55} { set col 2 }\n")
        f.write("        if {$cy >= -0.01 && $cy <= 0.74676} { set row 0 } elseif {$cy > 0.74676 && $cy <= 1.49352} { set row 1 } elseif {$cy > 1.49352 && $cy <= 2.25} { set row 2 }\n")
        
        f.write("        if {$col != -1 && $row != -1} {\n")
        f.write("            set panel_id [expr {$col * 3 + $row + 1}]\n")
        f.write("            if {$is_upper == 1} {\n")
        f.write("                set panel_id [expr {$panel_id + 9}]\n")
        f.write("                lappend list_Upper_Panel_$panel_id $elem_id\n")
        f.write("            } else {\n")
        f.write("                lappend list_Lower_Panel_$panel_id $elem_id\n")
        f.write("            }\n")
        f.write("        } else {\n")
        f.write("            lappend list_Spars_and_Ribs $elem_id\n")
        f.write("        }\n")
        f.write("    } else {\n")
        f.write("        lappend list_Spars_and_Ribs $elem_id\n")
        f.write("    }\n")
        f.write("}\n\n")
        
        for zone in zones: 
            f.write(f"if {{[llength $list_{zone['component']}] > 0}} {{\n    eval *createmark elems 1 $list_{zone['component']}\n    catch {{ *movemark elems 1 \"{zone['component']}\" }}\n}}\n")
        f.write("if {[llength $list_Spars_and_Ribs] > 0} {\n    eval *createmark elems 1 $list_Spars_and_Ribs\n    catch { *movemark elems 1 \"Spars_and_Ribs\" }\n}\n\n")

        # 4. EXPORT PANEL MAP FILE AND MATERIAL DEFINITIONS
        map_path_tcl = "C:/Users/Hugo/OneDrive/Escritorio/sim_part1/panel_mapping.txt"
        f.write(f"set fp [open \"{map_path_tcl}\" w]\n")
        for zone in zones: 
            f.write(f"puts $fp \"{zone['component']}: $list_{zone['component']}\"\n")
        f.write("close $fp\n\n")

        comp = materials["Composite"]
        f.write(f'if {{[hm_entityinfo exist mats "{comp["name"]}" -byname] == 0}} {{ *createentity mats cardimage="MAT8" name="{comp["name"]}" }}\n')
        f.write(f'*createmark mats 1 "{comp["name"]}"\n')
        f.write("set mat_c_id [lindex [hm_getmark mats 1] 0]\n")
        f.write(f"*setvalue mats id=$mat_c_id STATUS=2 196={comp['E1']} 197={comp['E2']} 198={comp['NU12']} 199={comp['G12']} 202={comp['RHO']}\n\n")

        alu = materials["Aluminum"]
        f.write(f'if {{[hm_entityinfo exist mats "{alu["name"]}" -byname] == 0}} {{ *createentity mats cardimage="MAT1" name="{alu["name"]}" }}\n')
        f.write(f'*createmark mats 1 "{alu["name"]}"\n')
        f.write("set mat_a_id [lindex [hm_getmark mats 1] 0]\n")
        f.write(f"*setvalue mats id=$mat_a_id STATUS=2 1={alu['E']} 3={alu['NU']} 4={alu['RHO']}\n\n")

        # 5. ASSIGN PCOMP PROPERTIES TO SKIN COVER PANELS
        created_pcomps = set()
        for zone in zones:
            prop_name = zone['name']
            if prop_name not in created_pcomps:
                f.write(f"*createentity props cardimage=\"PCOMP\" name=\"{prop_name}\"\n")
                f.write(f"*createmark props 1 \"{prop_name}\"\nset prop_id [lindex [hm_getmark props 1] 0]\n")
                num_plies = len(zone['stack'])
                if num_plies > 0:
                    z0_offset = -0.5 * (num_plies * zone['thickness'])
                    mats_str = " ".join(["$mat_c_id"] * num_plies)
                    thicks_str = " ".join([str(zone['thickness'])] * num_plies)
                    angs_str = " ".join([str(a) for a in zone['stack']])
                    f.write(f"eval *setvalue props id=$prop_id STATUS=2 3027={num_plies} 3014={z0_offset} 3023={{mats {mats_str}}} 3024={{{thicks_str}}} 3025={{{angs_str}}}\n")
                created_pcomps.add(prop_name)
            f.write(f"*createmark comps 1 \"{zone['component']}\"\nset comp_id [lindex [hm_getmark comps 1] 0]\n")
            f.write("catch { *setvalue comps id=$comp_id propertyid=$prop_id }\n")
            f.write(f"*createmark elems 1 \"by comp\" \"{zone['component']}\"\nif {{[hm_marklength elems 1] > 0}} {{ catch {{ *setvalue elems mark=1 propertyid=$prop_id }} }}\n\n")

        # 6. ASSIGN PSHELL TO SPARS AND RIBS
        f.write("\n# ---> SPARS AND RIBS PROPERTY ASSIGNMENT <--- \n")
        prop_name_spars = "PSHELL_Spars_Ribs"
        target_spars_plies = spars_plies_override if spars_plies_override is not None else 44
        alu_thickness = target_spars_plies * PLY_THICKNESS_SI 
        f.write(f"*createentity props cardimage=\"PSHELL\" name=\"{prop_name_spars}\"\n")
        f.write(f"*createmark props 1 \"{prop_name_spars}\"\nset pshell_id [lindex [hm_getmark props 1] 0]\n")
        f.write("*setvalue props id=$pshell_id materialid=$mat_a_id\n")
        f.write(f"*setvalue props id=$pshell_id STATUS=2 95={alu_thickness}\n")
        f.write("*createmark comps 1 \"Spars_and_Ribs\"\nset comp_id [lindex [hm_getmark comps 1] 0]\n")
        f.write("catch { *setvalue comps id=$comp_id propertyid=$pshell_id }\n")
        f.write("*createmark elems 1 \"by comp\" \"Spars_and_Ribs\"\nif {[hm_marklength elems 1] > 0} { catch { *setvalue elems mark=1 propertyid=$pshell_id } }\n\n")

        # 7. LOADS AND BOUNDARY CONDITIONS
        f.write("# --- Identify Root (min_x) and Tip (max_x) Extremes ---\n")
        f.write("*createmark nodes 1 all\nset min_x 99999.0\nset max_x -99999.0\n")
        f.write("if {[hm_marklength nodes 1] > 0} {\n    foreach n [hm_getmark nodes 1] {\n")
        f.write("        set x [hm_getvalue nodes id=$n dataname=x]\n")
        f.write("        if {$x < $min_x} { set min_x $x }\n        if {$x > $max_x} { set max_x $x }\n    }\n}\n\n")

        # Fixed boundary conditions (SPC) at root
        f.write("# --- SPC (Root Clamp) ---\n")
        f.write("if {[hm_entityinfo exist loadcols \"SPC\" -byname] == 0} { *createentity loadcols name=\"SPC\" }\n")
        f.write("set root_x0 [expr {$min_x - 0.0127}]; set root_x1 [expr {$min_x + 0.0127}]\n")
        f.write("*createmark nodes 1 \"by box\" $root_x0 -999 -999 $root_x1 999 999 0 inside 0 1 0\n")
        f.write("if {[hm_marklength nodes 1] > 0} {\n    foreach n [hm_getmark nodes 1] {\n        *createmark nodes 1 $n\n        catch { *loadcreate 1 3 1 0 0 0 0 0 0 }\n    }\n}\n")
        f.write("catch { *createmark loads 1 all; *movemark loads 1 \"SPC\" }\n\n")

        # Point loads applied at tip mid-height (Z_mid = 0.1905 m)
        f.write("# --- Point Forces (Tip) Applied at Z_mid = 0.1905 m ---\n")
        f.write("if {[hm_entityinfo exist loadcols \"Forces\" -byname] == 0} { *createentity loadcols name=\"Forces\" }\n")
        f.write("*createmark loadcols 1 \"Forces\"\nset force_id [lindex [hm_getmark loadcols 1] 0]\n\n")
        
        point_loads = [ 
            {"y": 0.0,     "Fz_N": 20235.01}, 
            {"y": 0.74676, "Fz_N": 42239.02}, 
            {"y": 1.49352, "Fz_N": 42239.02}, 
            {"y": 2.24028, "Fz_N": 85467.03} 
        ]
        f.write("set tip_x0 [expr {$max_x - 0.0254}]; set tip_x1 [expr {$max_x + 0.0254}]\n")
        
        for load in point_loads:
            f.write(f"set y0 [expr {{{load['y']} - 0.0508}}]; set y1 [expr {{{load['y']} + 0.0508}}]\n")
            f.write("*createmark nodes 1 \"by box\" $tip_x0 $y0 -0.0254 $tip_x1 $y1 0.50 0 inside 0 1 0\n")
            f.write("set candidate_nodes [hm_getmark nodes 1]\nset best_node -1\nset min_dz 9999.0\n")
            f.write("if {[llength $candidate_nodes] > 0} {\n    foreach n $candidate_nodes {\n")
            f.write("        set z [hm_getvalue nodes id=$n dataname=z]\n")
            f.write(f"        set dz [expr {{abs($z - {Z_MID})}}]\n")  # Select node closest to Z_MID = 0.1905 m
            f.write("        if {$dz < $min_dz} { set min_dz $dz; set best_node $n }\n    }\n")
            f.write("    if {$best_node != -1} {\n        *createmark nodes 1 $best_node\n")
            f.write(f"        catch {{ *loadcreate 1 1 1 0.0 0.0 {load['Fz_N']} 0.0 0.0 0.0 }}\n    }}\n}}\n")
            
        f.write("catch { *createmark loads 1 \"by box\" $tip_x0 -10.0 -10.0 $tip_x1 100.0 30.0 0 inside 0 1 0; *movemark loads 1 \"Forces\" }\n\n")
        f.write("*createmark loadcols 1 \"SPC\"\nset spc_id [lindex [hm_getmark loadcols 1] 0]\n\n")

        # 8. SUBCASE: LINEAR STATIC
        f.write("# --- 8. SUBCASE: LINEAR STATIC ---\n")
        f.write("if {[hm_entityinfo exist loadsteps \"EVAL_STATIC\" -byname] == 0} { *createentity loadsteps name=\"EVAL_STATIC\" }\n")
        f.write("*createmark loadsteps 1 \"EVAL_STATIC\"\nset ls_static_id [lindex [hm_getmark loadsteps 1] 0]\n\n")
        f.write("*attributeupdateint loadsteps $ls_static_id 4143 1 1 0 1\n*attributeupdateentity loadsteps $ls_static_id 4145 1 1 0 loadcols $spc_id\n")
        f.write("*attributeupdateint loadsteps $ls_static_id 4709 1 1 0 1\n*attributeupdateentity loadsteps $ls_static_id 4147 1 1 0 loadcols $force_id\n")
        f.write("*setvalue loadsteps id=$ls_static_id STATUS=2 4059=1 4060=STATICS\n\n")

        # 9. SUBCASE: BUCKLING
        f.write("# --- 9. SUBCASE: BUCKLING ---\n")
        f.write("if {[hm_entityinfo exist loadstepinputs \"EIGRL_BUCK\" -byname] == 0} { *createentity loadstepinputs cardimage=EIGRL name=\"EIGRL_BUCK\" }\n")
        f.write("*createmark loadstepinputs 1 \"EIGRL_BUCK\"\nset eigrl_id [lindex [hm_getmark loadstepinputs 1] 0]\n")
        f.write("catch { *setvalue loadstepinputs id=$eigrl_id STATUS=1 298=2 }\n\n")

        f.write("if {[hm_entityinfo exist loadsteps \"EVAL_BUCKLING\" -byname] == 0} { *createentity loadsteps name=\"EVAL_BUCKLING\" }\n")
        f.write("*createmark loadsteps 1 \"EVAL_BUCKLING\"\nset ls_buck_id [lindex [hm_getmark loadsteps 1] 0]\n\n")
        f.write("*attributeupdateint loadsteps $ls_buck_id 4143 1 1 0 1\n*attributeupdateentity loadsteps $ls_buck_id 4145 1 1 0 loadcols $spc_id\n")
        f.write("*attributeupdateint loadsteps $ls_buck_id 4709 1 1 0 1\n*attributeupdateentity loadsteps $ls_buck_id 4711 1 1 0 loadsteps $ls_static_id\n")
        f.write("*attributeupdateint loadsteps $ls_buck_id 4186 1 1 0 1\n*attributeupdateentity loadsteps $ls_buck_id 4187 1 1 0 loadstepinputs $eigrl_id\n")
        f.write("*setvalue loadsteps id=$ls_buck_id STATUS=2 4059=1 4060=BUCKLING\n\n")

        # 10. CONTROL CARDS & OUTPUT REQUESTS
        f.write("# --- Control Card Configuration ---\n")
        f.write("catch { *createentity controlcards cardimage=GLOBAL_OUTPUT_REQUEST }\n")
        f.write("catch { *createentity controlcards cardimage=UNSUPPORTED_CONTROL_CARDS }\n\n")
        
        # Inject Nastran/OptiStruct parameters
        output_req = (
            "PARAM,PRTMSS,YES\\n"            # Enables mass breakdown per PCOMP/PSHELL/Component
            "ESE(H3D, PRINT, PERCENT) = ALL\\n" # Requests Strain Energy in H3D, .out, and percentage distribution
            "FORCE(PUNCH) = ALL"             # Requests force exports to PUNCH file
        )
        f.write(f"set output_req \"{output_req}\"\n")
        f.write("catch { *setvalue controlcards cardimage=UNSUPPORTED_CONTROL_CARDS text=$output_req }\n")
        f.write("hm_blockmessages 0\nputs \"\\n---> MODEL PROCESSED SUCCESSFULLY (Mass and ESE enabled). <---\\n\"\n")

    return zones

# =================================================================
# --- 3. POST-PROCESSING
# =================================================================
def evaluate_results_from_punch(zones):
    """
    Parses internal forces from the OptiStruct PUNCH (.pch) file and calculates 
    analytical buckling load factors for each upper panel zone.
    """
    PCH_FILE = r"C:\Users\Hugo\OneDrive\Escritorio\sim_part1\sim1.pch" 
    MAP_FILE = r"C:\Users\Hugo\OneDrive\Escritorio\sim_part1\panel_mapping.txt" 
    
    if not os.path.exists(PCH_FILE) or not os.path.exists(MAP_FILE):
        print("\n[INFO] Execute TCL in HyperMesh and run OptiStruct before extracting results.")
        return

    print("\n[✓] Reading element mapping file...")
    panel_elements = {}
    with open(MAP_FILE, "r") as f:
        for line in f:
            if ":" in line:
                comp, elems = line.split(":")
                panel_elements[comp.strip()] = [int(e) for e in elems.strip().split() if e.strip().isdigit()]

    print("[✓] Extracting in-plane forces from PUNCH file (.pch)...")
    element_forces = {}
    current_label = ""
    reading_forces = False
    
    with open(PCH_FILE, "r") as f:
        for line in f:
            if line.startswith("$"):
                header = line.upper()
                if "LABEL" in header: 
                    current_label = header
                if "ELEMENT FORCES" in header:
                    if "STATIC" in current_label or "SUBCASE 1" in current_label: 
                        reading_forces = True
                    else: 
                        reading_forces = False 
                continue
            
            if reading_forces:
                if line.startswith("-CONT-"): 
                    continue
                tokens = line.strip().split()
                if len(tokens) >= 4 and tokens[0].isdigit():
                    elem_id = int(tokens[0])
                    try:
                        fx = float(tokens[1])
                        fy = float(tokens[2])
                        fxy = float(tokens[3])
                        element_forces[elem_id] = (fx, fy, fxy)
                    except ValueError: 
                        pass

    print("\n" + "="*60 + "\n   ANALYTICAL BUCKLING RESULTS IN SI UNITS (Upper Panels)\n" + "="*60)
    for zone in zones:
        comp_name = zone['component']
        elems = panel_elements.get(comp_name, [])
        if not elems: 
            continue

        sum_nx = sum_ny = sum_nxy = 0.0
        count = 0
        for eid in elems:
            if eid in element_forces:
                fx, fy, fxy = element_forces[eid]
                sum_nx += abs(fx)
                sum_ny += abs(fy)
                sum_nxy += fxy
                count += 1
                
        if count == 0: 
            continue
        Nx = sum_nx / count
        Ny = sum_ny / count
        Nxy = sum_nxy / count
        
        d11, d22, d12, d66 = calculate_D_matrix(
            E1=127.55e9, E2=13.03e9, nu12=0.3, G12=6.412e9, 
            stack=zone['stack'], ply_thickness=zone['thickness']
        )
        
        lam_c = calculate_analytical_buckling(d11, d22, d12, d66, 1.1811, 0.74676, Nx, Ny, Nxy)
        print(f"[{comp_name}] | Nx={Nx:9.1f} N/m | Ny={Ny:9.1f} N/m | Nxy={Nxy:9.1f} N/m | Buckling Lambda: {lam_c:.4f}")
    print("="*60 + "\n")


def apply_fem_parameters(fem_path="sim1.fem"):
    """
    Modifies an OptiStruct/Nastran input deck (.fem) to inject PARAM, PRTMSS, YES 
    and request Strain Energy percentage output.
    """
    if not os.path.exists(fem_path):
        print(f"[!] File not found: {fem_path}")
        return

    with open(fem_path, "r", encoding="utf-8") as f:
        content = f.readlines()

    new_content = []
    for line in content:
        # 1. Modify ESE output request to include PRINT and PERCENT
        if line.strip().startswith("ESE =") or line.strip().startswith("ESE("):
            new_content.append("ESE(PRINT, PERCENT) = ALL\n")
        # 2. Insert PARAM, PRTMSS, YES directly after BEGIN BULK
        elif line.strip() == "BEGIN BULK":
            new_content.append(line)
            new_content.append("PARAM, PRTMSS, YES\n")
        else:
            new_content.append(line)

    with open(fem_path, "w", encoding="utf-8") as f:
        f.writelines(new_content)
        
    print("[✓] .fem deck updated with PARAM, PRTMSS, YES and ESE(PRINT, PERCENT).")


if __name__ == "__main__":    
    lower_repetitions = {0: 26, 1: 21, 2: 15}
    upper_repetitions = {0: 26, 1: 21, 2: 15}
    
    generated_zones = generate_wing_box_model(
        dd_angles=[9.0, 80.0],
        lower_reps=lower_repetitions,
        upper_reps=upper_repetitions
    )

    evaluate_results_from_punch(generated_zones)
    apply_fem_parameters("sim1.fem")