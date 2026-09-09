
import os
import math

# =============================================================================
# 1. AUXILIARY LAYUP GENERATION FUNCTIONS
# =============================================================================
def generate_dd_stack(total_plies, phi, psi):
    """
    Generates a Double-Double [phi, -psi, -phi, psi] stacking sequence.
    Calculates repetition count to match the requested ply count exactly.
    
    Args:
        total_plies (int): Total number of plies in the laminate.
        phi (float): Primary ply angle in degrees.
        psi (float): Secondary ply angle in degrees.
        
    Returns:
        list: Complete list of orientation angles.
    """
    if total_plies <= 0:
        return []
        
    block = [phi, -psi, -phi, psi]
    reps = total_plies // 4
    remainder = total_plies % 4
    
    # Repeat building block and append remainder plies if not divisible by 4
    stack = (block * reps) + block[:remainder]
    return stack


# =============================================================================
# 2. ANALYTICAL BUCKLING SOLVER (Liu 2000 / Seresta 2007)
# =============================================================================
def calculate_D_matrix(E1, E2, nu12, G12, stack, ply_thickness=0.005):
    """
    Calculates the flexural stiffness matrix [D] based on CLPT.
    
    Args:
        E1, E2, nu12, G12 (float): Elastic material constants.
        stack (list): Laminate stacking sequence of ply angles.
        ply_thickness (float): Individual ply thickness.
        
    Returns:
        tuple: Flexural stiffness components (D11, D22, D12, D66).
    """
    nu21 = (E2 / E1) * nu12
    Q11 = E1 / (1 - nu12 * nu21)
    Q22 = E2 / (1 - nu12 * nu21)
    Q12 = (nu12 * E2) / (1 - nu12 * nu21)
    Q66 = G12
    
    D11 = D22 = D12 = D66 = 0.0
    total_thickness = len(stack) * ply_thickness
    z_current = -total_thickness / 2.0
    
    for angle in stack:
        theta = math.radians(angle)
        m, n = math.cos(theta), math.sin(theta)
        m2, n2, m4, n4, m2n2 = m**2, n**2, m**4, n**4, (m**2)*(n**2)
        
        # Transformed reduced stiffnesses
        Q11_bar = Q11*m4 + 2*(Q12 + 2*Q66)*m2n2 + Q22*n4
        Q22_bar = Q11*n4 + 2*(Q12 + 2*Q66)*m2n2 + Q22*m4
        Q12_bar = (Q11 + Q22 - 4*Q66)*m2n2 + Q12*(m4 + n4)
        Q66_bar = (Q11 + Q22 - 2*Q12 - 2*Q66)*m2n2 + Q66*(m4 + n4)
        
        z_k_minus_1 = z_current
        z_k = z_current + ply_thickness
        z_factor = (z_k**3 - z_k_minus_1**3) / 3.0
        
        D11 += Q11_bar * z_factor
        D22 += Q22_bar * z_factor
        D12 += Q12_bar * z_factor
        D66 += Q66_bar * z_factor
        z_current = z_k
        
    return D11, D22, D12, D66


def calculate_analytical_buckling(D11, D22, D12, D66, a, b, Nx, Ny, Nxy):
    """
    Calculates combined critical buckling load factor (lambda_c) for a simply supported plate.
    
    Args:
        D11, D22, D12, D66 (float): Flexural stiffness terms.
        a, b (float): Plate length (X) and width (Y).
        Nx, Ny, Nxy (float): Applied membrane force intensities per unit length.
        
    Returns:
        float: Combined critical buckling factor lambda_c.
    """
    lambda_n_min = float('inf')
    
    # Biaxial compression interaction loop (Nx, Ny)
    for m in range(1, 30):
        for n in range(1, 15):
            num = (math.pi**2) * (D11 * (m/a)**4 + 2*(D12 + 2*D66) * ((m/a)**2) * ((n/b)**2) + D22 * (n/b)**4)
            den = ((m/a)**2) * Nx + ((n/b)**2) * Ny
            if den > 0:
                lambda_n = num / den
                if lambda_n < lambda_n_min: 
                    lambda_n_min = lambda_n

    if lambda_n_min == float('inf'): 
        lambda_n_min = 1e6 

    # In-plane shear buckling interaction (Nxy)
    Nxy_abs = abs(Nxy)
    if Nxy_abs < 1e-6:
        lambda_s = 1e6
    else:
        gamma = math.sqrt(D11 * D22) / (D12 + 2*D66)
        beta_1 = 10 
        if gamma >= 1.0: 
            lambda_s = (4 * beta_1 * (D11 * D22**3)**0.25) / (b**2 * Nxy_abs)
        else:            
            lambda_s = (4 * beta_1 * math.sqrt(D22 * (D12 + 2*D66))) / (b**2 * Nxy_abs)

    A = 1.0 / (lambda_s**2)
    B = 1.0 / lambda_n_min
    
    if A > 1e-12:
        # Positive root of quadratic formula for combined shear-compression buckling
        lambda_c = (-B + math.sqrt(B**2 + 4*A)) / (2*A)
    else:
        # Pure biaxial compression scenario
        lambda_c = lambda_n_min
        
    return lambda_c


# =============================================================================
# 3. HYPERMESH MODEL GENERATOR (PRE-PROCESSOR)
# =============================================================================
def generate_wing_box_model(
    stacking_type="traditional", 
    dd_angles=(30, 60), 
    dd_plies_override=None,
    spars_material="composite",        
    spars_stacking_type="traditional", 
    spars_dd_angles=(45, 45),
    spars_plies_override=None
):
    """
    Generates a TCL automation script for HyperMesh pre-processing.
    
    Configures skin panels (Upper/Lower), internal spars/ribs, materials (MAT8/MAT1),
    property cards (PCOMP/PSHELL), boundary conditions, point loads, and static subcases.
    """
    # Base model database path
    MODEL_PATH = r"c:\Users\Hugo\OneDrive\Escritorio\v2_paper\v2_mesh.hm" 
    OUTPUT_FILE = "v2.tcl"
    DEFAULT_MESH_SIZE = 2.5
    MESH_TYPE = 1  # Quads
    
    # Composite Material Properties (T300/5208 CFRP)
    E1_comp = 18.5e6; E2_comp = 1.89e6; G12_comp = 0.93e6; NU12 = 0.3
    
    materials = {
        "Composite": { "id": 1, "type": "MAT8", "name": "T300_5208", 
                       "E1": E1_comp, "E2": E2_comp, "NU12": NU12, "G12": G12_comp, "RHO": 0.057 },
        "Aluminum":  { "id": 2, "type": "MAT1", "name": "Al_2024_T3", "E": 10.0e6, "NU": 0.3, "RHO": 0.1 }
    }

    def exp(plies, reps): return plies * reps

    # Explicit baseline layup sequences for upper skin panels (10 to 18)
    explicit_upper_half_stacks = [
        exp([45, -45], 12) + exp([90], 4) + exp([0], 2) + exp([90, 90, 0, 0], 5) + exp([0, 0, 90, 90, 90, 90], 2) + exp([0], 4) + exp([90], 2) + exp([0], 2),
        exp([45, -45], 15) + exp([0], 2) + [45, -45] + exp([0], 4) + exp([90, 90, 90, 90, 0, 0], 3) + exp([90], 2),
        exp([45, -45], 12) + exp([0, 0, 90, 90, 90, 90], 2) + exp([0], 2) + exp([90], 2),
        exp([45, -45], 12) + exp([0], 2) + [45, -45] + exp([0], 4) + exp([90, 90, 90, 90, 0, 0], 2) + exp([90], 2),
        exp([45, -45], 16) + exp([0], 2) + [45, -45] + exp([0], 2) + exp([0, 0, 90, 90, 90, 90], 4) + exp([0], 2) + exp([90], 2),
        exp([45, -45], 11) + exp([90], 2) + exp([45, -45], 2) + exp([90], 2) + [45, -45] + exp([90], 2) + exp([90, 90, 0, 0], 3) + exp([0, 0, 90, 90, 0, 0], 6),
        exp([45, -45], 16) + exp([90, 90, 45, -45], 2) + exp([90, 90, 90, 90, 0, 0], 3) + exp([90, 90, 0, 0], 2) + exp([0, 0, 90, 90, 0, 0], 3),
        exp([45, -45], 20) + exp([0], 2) + [45, -45] + exp([0], 4) + exp([90, 90, 90, 90, 0, 0], 3) + exp([90], 2),
        exp([45, -45], 15) + exp([0], 4) + exp([90, 90, 0, 0], 2) + exp([90], 4) + exp([0], 2) + exp([90], 2)
    ]
    plies_lower = [(6, 1, 0), (4, 0, 0), (1, 0, 0), (1, 2, 0), (4, 2, 1), (7, 1, 0), (11, 1, 1), (7, 1, 1), (2, 2, 0)]

    zones = []
    x_coords = [0.0, 46.5, 93.0]; y_coords = [0.0, 29.4, 58.8]
    panel_grid_map = {
        1: (0, 0), 2: (1, 0), 3: (2, 0), 6: (0, 1), 5: (1, 1), 4: (2, 1), 7: (0, 2), 8: (1, 2), 9: (2, 2),
        10: (0, 0), 11: (1, 0), 12: (2, 0), 15: (0, 1), 14: (1, 1), 13: (2, 1), 16: (0, 2), 17: (1, 2), 18: (2, 2)
    }

    def get_dd_ply_count(p_id, auto_count):
        if dd_plies_override is None: return auto_count
        elif isinstance(dd_plies_override, int): return dd_plies_override
        elif isinstance(dd_plies_override, dict): return dd_plies_override.get(p_id, auto_count)
        return auto_count

    # --- LOWER SKIN PANELS (Panels 1 to 9) ---
    for panel_id in range(1, 10):
        col, row = panel_grid_map[panel_id]
        n0, n45, n90 = plies_lower[panel_id - 1]
        half = []; [half.extend([45, -45]) for _ in range(n45)]
        half.extend([0, 0] * n0); half.extend([90, 90] * n90)
        
        if stacking_type == "double_double":
            auto_plies = len(half) * 2
            target_plies = get_dd_ply_count(panel_id, auto_plies)
            final_stack = generate_dd_stack(target_plies, dd_angles[0], dd_angles[1])
        else:
            final_stack = half + half[::-1]
            
        zones.append({
            "name": f"PCOMP_Low_P{panel_id}", "component": f"Lower_Panel_{panel_id}", "stack": final_stack, 
            "thickness": 0.005, "x_min": x_coords[col], "x_max": x_coords[col] + 46.5, "y_min": y_coords[row], "y_max": y_coords[row] + 29.4
        })

    # --- UPPER SKIN PANELS (Panels 10 to 18) ---
    for panel_id in range(10, 19):
        col, row = panel_grid_map[panel_id]
        half = explicit_upper_half_stacks[panel_id - 10]
        
        if stacking_type == "double_double":
            auto_plies = len(half) * 2
            target_plies = get_dd_ply_count(panel_id, auto_plies)
            final_stack = generate_dd_stack(target_plies, dd_angles[0], dd_angles[1])
        else:
            final_stack = half + half[::-1]
        
        zones.append({
            "name": f"PCOMP_Up_P{panel_id}", "component": f"Upper_Panel_{panel_id}", "stack": final_stack, 
            "thickness": 0.005, "x_min": x_coords[col], "x_max": x_coords[col] + 46.5, "y_min": y_coords[row], "y_max": y_coords[row] + 29.4
        })

    print(f"[INFO] Generating TCL automation script... Skins: {stacking_type.upper()} | Spars: {spars_material.upper()}")
    
    with open(OUTPUT_FILE, "w", encoding='utf-8') as f:
        f.write("# HyperMesh Model Automation Engine - Generated via Python\n")
        f.write("hm_answernext yes\n")
        f.write("hm_blockmessages 1\n\n")

        # Step 1: Global Cleanup
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 1: Model Database Cleanup\n")
        f.write("# ---------------------------------------------------------\n")
        f.write("catch { *clearmark loads 1; *createmark loads 1 all; *deletemark loads 1 }\n")
        f.write("catch { *clearmark elems 1; *createmark elems 1 all; *deletemark elems 1 }\n")
        f.write("catch { *clearmark nodes 1; *createmark nodes 1 all; *deletemark nodes 1 }\n")
        f.write("catch { *clearmark props 1; *createmark props 1 all; *deletemark props 1 }\n")
        f.write("catch { *clearmark mats 1; *createmark mats 1 all; *deletemark mats 1 }\n\n")

        # Step 2: Surface Remeshing
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 2: Surface Remeshing & Node Equivalence\n")
        f.write("# ---------------------------------------------------------\n")
        f.write("*createmark surfs 1 all\n")
        f.write("catch { *edgeequivalence surfs 1 0.1 }\n")
        f.write(f"catch {{ *defaultremeshsurf 1 {DEFAULT_MESH_SIZE} {MESH_TYPE} 2 2 1 1 1 1 0 0 0 0 }}\n")
        f.write("*createmark nodes 1 all\n")
        f.write("if {[hm_marklength nodes 1] > 0} {\n")
        f.write("    catch { *equivalence nodes 1 1.0 1 0 0 }\n")
        f.write("}\n\n")

        # Step 3: Spatial Filtering & Element Routing
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 3: Spatial Sorting of Elements by Centroid (Z-Axis Filtering)\n")
        f.write("# ---------------------------------------------------------\n")
        for zone in zones: 
            f.write(f"if {{[hm_entityinfo exist comps \"{zone['component']}\" -byname] == 0}} {{\n")
            f.write(f"    *createentity comps name=\"{zone['component']}\"\n")
            f.write("}\n")
            
        f.write("if {[hm_entityinfo exist comps \"Spars_and_Ribs\" -byname] == 0} {\n")
        f.write("    *createentity comps name=\"Spars_and_Ribs\"\n")
        f.write("}\n\n")
        
        f.write("*createmark elems 1 all\n")
        f.write("set all_elems [hm_getmark elems 1]\n")
        for zone in zones: 
            f.write(f"set list_{zone['component']} \"\"\n")
        f.write("set list_Spars_and_Ribs \"\"\n\n")
        
        f.write("foreach elem_id $all_elems {\n")
        f.write("    set n [hm_getvalue elems id=$elem_id dataname=nodes]\n")
        f.write("    set cx [expr {([hm_getvalue nodes id=[lindex $n 0] dataname=x] + [hm_getvalue nodes id=[lindex $n 1] dataname=x] + [hm_getvalue nodes id=[lindex $n 2] dataname=x]) / 3.0}]\n")
        f.write("    set cy [expr {([hm_getvalue nodes id=[lindex $n 0] dataname=y] + [hm_getvalue nodes id=[lindex $n 1] dataname=y] + [hm_getvalue nodes id=[lindex $n 2] dataname=y]) / 3.0}]\n")
        f.write("    set cz [expr {([hm_getvalue nodes id=[lindex $n 0] dataname=z] + [hm_getvalue nodes id=[lindex $n 1] dataname=z] + [hm_getvalue nodes id=[lindex $n 2] dataname=z]) / 3.0}]\n")
        f.write("    set assigned 0\n")

        for i, zone in enumerate(zones):
            z_ideal = 15.0 if "Up" in zone['name'] else 0.0
            cond = f"$cx >= {zone['x_min'] - 0.1} && $cx <= {zone['x_max'] + 0.1} && $cy >= {zone['y_min'] - 0.1} && $cy <= {zone['y_max'] + 0.1} && $cz >= {z_ideal - 0.1} && $cz <= {z_ideal + 0.1}"
            
            if i == 0: f.write(f"    if {{{cond}}} {{\n")
            else:      f.write(f"    }} elseif {{{cond}}} {{\n")
                
            f.write(f"        lappend list_{zone['component']} $elem_id\n")
            f.write(f"        set assigned 1\n")
            
        f.write("    }\n")
        f.write("    if {$assigned == 0} {\n")
        f.write("        lappend list_Spars_and_Ribs $elem_id\n")
        f.write("    }\n")
        f.write("}\n\n")
        
        for zone in zones: 
            f.write(f"if {{[llength $list_{zone['component']}] > 0}} {{\n")
            f.write(f"    eval *createmark elems 1 $list_{zone['component']}\n")
            f.write(f"    catch {{ *movemark elems 1 \"{zone['component']}\" }}\n")
            f.write("}\n")
            
        f.write("if {[llength $list_Spars_and_Ribs] > 0} {\n")
        f.write("    eval *createmark elems 1 $list_Spars_and_Ribs\n")
        f.write("    catch { *movemark elems 1 \"Spars_and_Ribs\" }\n")
        f.write("}\n\n")

        # Export panel element mapping file for post-processing
        map_path_tcl = "C:/Users/Hugo/OneDrive/Escritorio/v2_paper/panel_mapping.txt"
        f.write(f"set fp [open \"{map_path_tcl}\" w]\n")
        for zone in zones:
            f.write(f"puts $fp \"{zone['component']}: $list_{zone['component']}\"\n")
        f.write("close $fp\n\n")

        # Step 4: Material Definition
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 4: Create Materials (MAT8 Composite & MAT1 Isotropic)\n")
        f.write("# ---------------------------------------------------------\n")
        comp = materials["Composite"]
        f.write(f'*createentity mats cardimage="MAT8" name="{comp["name"]}"\n')
        f.write(f'*createmark mats 1 "{comp["name"]}"\n')
        f.write("set m_id [lindex [hm_getmark mats 1] 0]\n")
        f.write(f"*setvalue mats id=$m_id STATUS=2 196={comp['E1']} 197={comp['E2']} 198={comp['NU12']} 199={comp['G12']} 202={comp['RHO']}\n\n")

        alu = materials["Aluminum"]
        f.write(f'*createentity mats cardimage="MAT1" name="{alu["name"]}"\n')
        f.write(f'*createmark mats 1 "{alu["name"]}"\n')
        f.write("set m_alu_id [lindex [hm_getmark mats 1] 0]\n")
        f.write(f"*setvalue mats id=$m_alu_id STATUS=2 1={alu['E']} 3={alu['NU']} 4={alu['RHO']}\n\n")

        # Step 5: Skin Panel Property Cards (PCOMP)
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 5: Create PCOMP Property Cards for Skin Panels\n")
        f.write("# ---------------------------------------------------------\n")
        for zone in zones:
            f.write(f"*createentity props cardimage=\"PCOMP\" name=\"{zone['name']}\"\n")
            f.write(f"*createmark props 1 \"{zone['name']}\"\n")
            f.write("set p_id [lindex [hm_getmark props 1] 0]\n")
            num_plies = len(zone['stack'])
            
            if num_plies > 0:
                total_thickness = num_plies * zone['thickness']
                z0_offset = -0.5 * total_thickness
                
                mats_s = " ".join(["$m_id"] * num_plies)
                thicks_s = " ".join([str(zone['thickness'])] * num_plies)
                angs_s = " ".join([str(a) for a in zone['stack']])
                
                f.write(f"*setvalue props id=$p_id STATUS=2 3027={num_plies} 3014={z0_offset} 3021=0.0 3023={{mats {mats_s}}} 3024={{{thicks_s}}} 3025={{{angs_s}}}\n")
                
                # Bind property card to component
                f.write(f"*createmark comps 1 \"{zone['component']}\"\n")
                f.write(f"catch {{ *propertyupdate comps 1 \"{zone['name']}\" }}\n\n")

        # Step 6: Spars & Ribs Property Setup
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 6: Spars & Ribs Property Configuration (PCOMP vs PSHELL)\n")
        f.write("# ---------------------------------------------------------\n")
        
        target_spars_plies = spars_plies_override if spars_plies_override is not None else 44
        prop_name_spars = ""
        
        if spars_material.lower() == "aluminum":
            prop_name_spars = "PSHELL_Spars_Ribs"
            f.write(f"*createentity props cardimage=\"PSHELL\" name=\"{prop_name_spars}\"\n")
            f.write(f"*createmark props 1 \"{prop_name_spars}\"\n")
            f.write("set p_id [lindex [hm_getmark props 1] 0]\n")
            
            alu_thickness = target_spars_plies * 0.005 
            f.write("*setvalue props id=$p_id materialid=$m_alu_id\n")
            f.write(f"*setvalue props id=$p_id STATUS=2 95={alu_thickness}\n")
            
        elif spars_material.lower() == "composite":
            prop_name_spars = "PCOMP_Spars_Ribs"
            f.write(f"*createentity props cardimage=\"PCOMP\" name=\"{prop_name_spars}\"\n")
            f.write(f"*createmark props 1 \"{prop_name_spars}\"\n")
            f.write("set p_id [lindex [hm_getmark props 1] 0]\n")
            
            if spars_stacking_type == "double_double":
                s_stack = generate_dd_stack(target_spars_plies, spars_dd_angles[0], spars_dd_angles[1])
            else:
                reps = target_spars_plies // 2
                remainder = target_spars_plies % 2
                s_stack = [45, -45] * reps
                if remainder != 0:
                    s_stack.append(45) 
                    
            total_spars_thickness = len(s_stack) * 0.005
            z0_spars_offset = -0.5 * total_spars_thickness
                
            m_s = " ".join(["$m_id"] * len(s_stack))
            t_s = " ".join(["0.005"] * len(s_stack))
            a_s = " ".join([str(a) for a in s_stack])
            
            f.write(f"*setvalue props id=$p_id STATUS=2 3027={len(s_stack)} 3014={z0_spars_offset} 3021=0.0 3023={{mats {m_s}}} 3024={{{t_s}}} 3025={{{a_s}}}\n")
            
        f.write("*createmark comps 1 \"Spars_and_Ribs\"\n")
        f.write(f"catch {{ *propertyupdate comps 1 \"{prop_name_spars}\" }}\n\n")

        # Step 7: Boundary Conditions, Loads & Subcase Card Setup
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 7: Apply Single Point Constraints (SPC) & Point Forces\n")
        f.write("# ---------------------------------------------------------\n")
        f.write("if {[hm_entityinfo exist loadcols \"SPC\" -byname] == 0} {\n")
        f.write("    *createentity loadcols name=\"SPC\"\n")
        f.write("}\n")
        f.write("*createmark nodes 1 \"by box\" -0.5 -999 -999 0.5 999 999 0 inside 0 1 0\n")
        f.write("if {[hm_marklength nodes 1] > 0} {\n")
        f.write("    foreach n [hm_getmark nodes 1] {\n")
        f.write("        *createmark nodes 1 $n\n")
        f.write("        catch { *loadcreate 1 3 1 0 0 0 0 0 0 }\n")
        f.write("    }\n")
        f.write("}\n")
        f.write("catch { *createmark loads 1 all; *movemark loads 1 \"SPC\" }\n\n")

        f.write("if {[hm_entityinfo exist loadcols \"Forces\" -byname] == 0} {\n")
        f.write("    *createentity loadcols name=\"Forces\"\n")
        f.write("}\n")
        f.write("*createmark loadcols 1 \"Forces\"\n")
        f.write("set force_id [lindex [hm_getmark loadcols 1] 0]\n\n")
        
        point_loads = [
            {"y": 0.0, "Fz_N": 20235.01}, 
            {"y": 29.4, "Fz_N": 42239.02}, 
            {"y": 58.8, "Fz_N": 42239.02}, 
            {"y": 88.2, "Fz_N": 85467.03}
        ]
        
        for load in point_loads:
            f.write(f"*createmark nodes 1 \"by box\" 139.0 {load['y']-2.0} -1.0 140.0 {load['y']+2.0} 16.0 0 inside 0 1 0\n")
            f.write("set candidate_nodes [hm_getmark nodes 1]\n")
            f.write("set best_node -1\n")
            f.write("set min_dz 9999.0\n")
            f.write("if {[llength $candidate_nodes] > 0} {\n")
            f.write("    foreach n $candidate_nodes {\n")
            f.write("        set z [hm_getvalue nodes id=$n dataname=z]\n")
            f.write("        set dz [expr {abs($z - 7.5)}]\n")
            f.write("        if {$dz < $min_dz} {\n")
            f.write("            set min_dz $dz\n")
            f.write("            set best_node $n\n")
            f.write("        }\n")
            f.write("    }\n")
            f.write("    if {$best_node != -1} {\n")
            f.write("        *createmark nodes 1 $best_node\n")
            f.write(f"        catch {{ *loadcreate 1 1 1 0.0 0.0 {load['Fz_N']} 0.0 0.0 0.0 }}\n")
            f.write("    }\n")
            f.write("}\n\n")
            
        f.write("*createmark loads 1 \"by box\" 130.0 -10.0 -10.0 150.0 100.0 30.0 0 inside 0 1 0\n")
        f.write("catch { *movemark loads 1 \"Forces\" }\n\n")

        f.write("if {[hm_entityinfo exist controlcards \"GLOBAL_OUTPUT_REQUEST\" -byname] == 0} {\n")
        f.write("    *createentity controlcards cardimage=\"GLOBAL_OUTPUT_REQUEST\" name=\"GLOBAL_OUTPUT_REQUEST\"\n")
        f.write("}\n")
        f.write("*setvalue controlcards name=\"GLOBAL_OUTPUT_REQUEST\" STATUS=1 107=1 108=3\n\n")
        
        f.write("if {[hm_entityinfo exist loadsteps \"EVAL_STATIC\" -byname] == 0} {\n")
        f.write("    *createentity loadsteps name=\"EVAL_STATIC\"\n")
        f.write("}\n")
        f.write("*createmark loadsteps 1 \"EVAL_STATIC\"\n")
        f.write("set ls_id [lindex [hm_getmark loadsteps 1] 0]\n\n")

        f.write("*createmark loadcols 1 \"SPC\"\n")
        f.write("set spc_id [lindex [hm_getmark loadcols 1] 0]\n")
        f.write("*createmark loadcols 1 \"Forces\"\n")
        f.write("set force_id [lindex [hm_getmark loadcols 1] 0]\n\n")

        f.write("*attributeupdateint loadsteps $ls_id 4143 1 1 0 1\n")
        f.write("*attributeupdateentity loadsteps $ls_id 4145 1 1 0 loadcols $spc_id\n")
        f.write("*attributeupdateint loadsteps $ls_id 4709 1 1 0 1\n")
        f.write("*attributeupdateentity loadsteps $ls_id 4147 1 1 0 loadcols $force_id\n")
        f.write("*setvalue loadsteps id=$ls_id STATUS=2 4059=1 4060=STATICS\n\n")
        
        f.write("hm_blockmessages 0\n")
        f.write("puts \"\\n===================================================\"\n")
        f.write("puts \"[SUCCESS] MODEL READY FOR OPTISTRUCT SOLVER LAUNCH\"\n")
        f.write("puts \"===================================================\\n\"\n")

    return zones


# =============================================================================
# 4. POST-PROCESSOR (PUNCH .pch Parser & Analytical Buckling Extraction)
# =============================================================================
def evaluate_results_from_punch(zones):
    """
    Parses membrane element forces (Nx, Ny, Nxy) from OptiStruct PUNCH (.pch) output files
    and evaluates analytical buckling factors across upper skin panels.
    """
    PCH_FILE = r"C:\Users\Hugo\OneDrive\Escritorio\v2_paper\v2.pch" 
    MAP_FILE = r"C:\Users\Hugo\OneDrive\Escritorio\v2_paper\panel_mapping.txt" 
    
    if not os.path.exists(PCH_FILE) or not os.path.exists(MAP_FILE):
        print("\n[INFO] Run the generated TCL in HyperMesh and execute OptiStruct before extracting results.")
        return

    print("\n[INFO] Reading panel element mapping file...")
    panel_elements = {}
    with open(MAP_FILE, "r") as f:
        for line in f:
            if ":" in line:
                comp, elems = line.split(":")
                panel_elements[comp.strip()] = [int(e) for e in elems.strip().split() if e.strip().isdigit()]

    print("[INFO] Extracting in-plane membrane forces from PUNCH file (.pch)...")
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
                    if "STATIC" in current_label: 
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

    print(f"[DEBUG] Extracted membrane force tensors for {len(element_forces)} elements.")

    print("\n" + "="*70)
    print("   ANALYTICAL BUCKLING MARGIN RESULTS (Upper Compression Panels)")
    print("="*70)
    
    for zone in zones:
        comp_name = zone['component']
        if "Upper" not in comp_name: 
            continue
        
        elems = panel_elements.get(comp_name, [])
        if not elems: 
            continue

        sum_nx = sum_ny = sum_nxy = 0.0
        count = 0
        for eid in elems:
            if eid in element_forces:
                fx, fy, fxy = element_forces[eid]
                sum_nx += abs(fy)  
                sum_ny += abs(fx)
                sum_nxy += abs(fxy)
                count += 1
                
        if count == 0: 
            continue
        
        Nx = sum_nx / count
        Ny = sum_ny / count
        Nxy = sum_nxy / count
        
        d11, d22, d12, d66 = calculate_D_matrix(18.5e6, 1.89e6, 0.3, 0.93e6, zone['stack'])
        lam_c = calculate_analytical_buckling(d11, d22, d12, d66, 46.5, 29.4, Nx, Ny, Nxy)
        
        print(f"[{comp_name}] | Nx={Nx:7.1f} | Ny={Ny:7.1f} | Nxy={Nxy:7.1f} | Buckling Factor (Lambda_c): {lam_c:.4f}")
    
    print("="*70 + "\n")


if __name__ == "__main__":
    
    # --- SKIN PANEL CONFIGURATION ---
    panel_stacking_type = "traditional"  # Options: "traditional" or "double_double"
    panel_dd_angles = (15, 75)            # Double-Double angle parameters (phi, psi)
    panel_plies_override = None          # Optional ply override count
    
    # --- SPARS & RIBS CONFIGURATION ---
    spars_mat = "aluminum"               # Options: "composite" or "aluminum"
    spars_stacking_type = "traditional"  # Options: "traditional" or "double_double"
    spars_dd_angles = (15, 45)
    spars_plies_override = 44
    
    zones_created = generate_wing_box_model(
        stacking_type=panel_stacking_type, 
        dd_angles=panel_dd_angles, 
        dd_plies_override=panel_plies_override,
        
        spars_material=spars_mat,
        spars_stacking_type=spars_stacking_type,
        spars_dd_angles=spars_dd_angles,
        spars_plies_override=spars_plies_override
    )
    
    evaluate_results_from_punch(zones_created)
