
import os
import subprocess

def generate_multi_patch_final():
    """
    Generates and executes 'v1.tcl' to automate a multi-patch composite FEA setup.
    """
    # -------------------------------------------------------------------------
    # 1. FILE PATH CONFIGURATION
    # -------------------------------------------------------------------------
    # Path to Altair HyperMesh executable (2025.1 Win64)
    HM_EXE = r"C:\Program Files\Altair\2025.1\hwdesktop\hm\bin\win64\hmopengl.exe"
    
    # Base geometry/mesh HyperMesh database
    MODEL_PATH = r"c:\Users\Hugo\OneDrive\Escritorio\v1_offset\v1.hm"
    
    # Target TCL automation script output
    OUTPUT_FILE = "v1.tcl"
    
    # -------------------------------------------------------------------------
    # 2. MATERIAL DEFINITIONS (MAT8 - Orthotropic 2D Shell Materials)
    # -------------------------------------------------------------------------
    # Units: Stiffness in MPa, Density in tonnes/mm^3
    materials = {
        1: {"name": "Carbon_1",    "E1": 140000, "E2": 10000, "NU12": 0.30, "G12": 5000, "RHO": 1.6e-9},
        2: {"name": "Carbon_2",    "E1": 180000, "E2": 12000, "NU12": 0.40, "G12": 7000, "RHO": 1.3e-9},
        3: {"name": "Glass_Fiber", "E1": 45000,  "E2": 10000, "NU12": 0.28, "G12": 4000, "RHO": 1.9e-9},
    }
    
    # -------------------------------------------------------------------------
    # 3. SPATIAL ZONE CONFIGURATION (100x50 mm Plate Partitioning)
    # -------------------------------------------------------------------------
    # Plate dimensions: 100 mm (X-axis) x 50 mm (Y-axis) divided into 4 quadrants
    zones = [
        # Quadrant 1: Bottom-Left (X: 0 to 50 mm | Y: 0 to 25 mm)
        {
            "name": "P_Patch_Q1", "component": "Zone_Q1", "material_id": 1,
            "stack": [25.0, -65.0, -25.0, 65.0], "thickness": 0.2, 
            "x_min": 0.0, "x_max": 50.0, "y_min": 0.0, "y_max": 25.0, "mesh_type": "CQUAD4"
        },
        # Quadrant 2: Bottom-Right (X: 50 to 100 mm | Y: 0 to 25 mm)
        {
            "name": "P_Patch_Q2", "component": "Zone_Q2", "material_id": 2,
            "stack": [15.0, -55.0, -15.0, 55.0], "thickness": 0.15, 
            "x_min": 50.0, "x_max": 100.0, "y_min": 0.0, "y_max": 25.0, "mesh_type": "CQUAD4"
        },
        # Quadrant 3: Top-Left (X: 0 to 50 mm | Y: 25 to 50 mm)
        {
            "name": "P_Patch_Q3", "component": "Zone_Q3", "material_id": 3,
            "stack": [35.0, -55.0, -35.0, 55.0], "thickness": 0.1, 
            "x_min": 0.0, "x_max": 50.0, "y_min": 25.0, "y_max": 50.0, "mesh_type": "CQUAD8"
        },
        # Quadrant 4: Top-Right (X: 50 to 100 mm | Y: 25 to 50 mm)
        {
            "name": "P_Patch_Q4", "component": "Zone_Q4", "material_id": 1,
            "stack": [45.0, -45.0, -45.0, 45.0], "thickness": 0.1, 
            "x_min": 50.0, "x_max": 100.0, "y_min": 25.0, "y_max": 50.0, "mesh_type": "CTRIA3"
        }
    ]

    # -------------------------------------------------------------------------
    # 4. BOUNDARY CONDITIONS (SPC) AND FORCE LOADS CONFIGURATION
    # -------------------------------------------------------------------------
    # Single Point Constraints (Fixed root at X = 0 mm, constrained DOFs 1-6)
    bcs = [
        {
            "name": "Fixed_Support_X0",
            "x_min": -0.1, "x_max": 0.1, 
            "y_min": -999.0, "y_max": 999.0, # Spans full Y height (0 to 50 mm)
            "z_min": -999.0, "z_max": 999.0,
            "dofs_values": "0.0 0.0 0.0 0.0 0.0 0.0" # Constrain UX, UY, UZ, RX, RY, RZ
        }
    ]

    # Concentrated vertical force load at tip (X = 100 mm, Fz = -50 N)
    loads = [
        {
            "name": "Bending_Force_Z_X100",
            "x_min": 99.9, "x_max": 100.1,
            "y_min": -999.0, "y_max": 999.0, # Spans full Y height (0 to 50 mm)
            "z_min": -999.0, "z_max": 999.0,
            "vector": "0.0 0.0 -50.0 0.0 0.0 0.0" # Fz = -50 N
        }
    ]
    
    # Format Windows path for TCL scripting engine
    tcl_model_path = MODEL_PATH.replace("\\", "/")
    print(f"[INFO] Generating multi-patch automation script: {OUTPUT_FILE}...")

    # -------------------------------------------------------------------------
    # 5. SCRIPT GENERATION ENGINE
    # -------------------------------------------------------------------------
    with open(OUTPUT_FILE, "w", encoding='utf-8') as f:
        f.write("# Dynamic Multi-Patch Composite FEA Setup - HyperMesh Automation\n")
        f.write("hm_answernext yes\n\n")
        
        # Load hypermesh database
        f.write("# ---------------------------------------------------------\n")
        f.write("# Load HyperMesh Model & Count Initial Mesh\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f'*readfile "{tcl_model_path}" 0\n\n')
        f.write("*createmark elems 1 all\n")
        f.write('set total_initial [hm_marklength elems 1]\n')
        f.write('puts "[INFO] Initial mesh element count: $total_initial"\n\n')

        # Section A: Component Creation
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section A: Create Target Zone Components\n")
        f.write("# ---------------------------------------------------------\n")
        for zone in zones:
            comp_name = zone["component"]
            f.write(f'if {{[hm_entityinfo exist comps "{comp_name}" -byname] == 0}} {{\n')
            f.write(f'    *createentity comps cardimage="" name="{comp_name}"\n')
            f.write(f'    puts "[INFO] Component {comp_name} created."\n')
            f.write(f'}} else {{\n')
            f.write(f'    puts "[INFO] Component {comp_name} already exists."\n')
            f.write(f'}}\n\n')

        # Section B: Dynamic Element Partitioning by 2D Centroid
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section B: Dynamic 2D Centroid Calculation & Element Sorting\n")
        f.write("# ---------------------------------------------------------\n")
        f.write("*createmark elems 1 all\n")
        f.write("set all_elems [hm_getmark elems 1]\n")
        f.write('puts "[INFO] Elements ready for spatial redistribution: [llength $all_elems]"\n\n')
        
        for i in range(len(zones)):
            f.write(f"set zone{i}_count 0\n")
        f.write("\n")
        
        # Loop through each element, compute mean nodal X/Y coordinates (centroid)
        f.write("foreach elem_id $all_elems {\n")
        f.write("    set nodes [hm_getvalue elems id=$elem_id dataname=nodes]\n")
        f.write("    set x_sum 0.0\n")
        f.write("    set y_sum 0.0\n")
        f.write("    set node_count 0\n")
        f.write("    foreach node_id $nodes {\n")
        f.write("        set x [hm_getvalue nodes id=$node_id dataname=x]\n")
        f.write("        set y [hm_getvalue nodes id=$node_id dataname=y]\n")
        f.write("        set x_sum [expr $x_sum + $x]\n")
        f.write("        set y_sum [expr $y_sum + $y]\n")
        f.write("        incr node_count\n")
        f.write("    }\n")
        f.write("    set centroid_x [expr $x_sum / $node_count]\n")
        f.write("    set centroid_y [expr $y_sum / $node_count]\n")
        f.write("    *createmark elems 1 $elem_id\n")
        
        # Spatial sorting into quadrant components based on bounding limits
        for i, zone in enumerate(zones):
            condicion = (f"$centroid_x >= {zone['x_min']} && $centroid_x <= {zone['x_max']} && "
                         f"$centroid_y >= {zone['y_min']} && $centroid_y <= {zone['y_max']}")
            if i == 0:
                f.write(f"    if {{{condicion}}} {{\n")
            else:
                f.write(f"    }} elseif {{{condicion}}} {{\n")
            f.write(f'        *movemark elems 1 "{zone["component"]}"\n')
            f.write(f"        incr zone{i}_count\n")
        
        f.write("    } else {\n")
        f.write('        *deletemark elems 1\n')  # Delete elements outside defined domains
        f.write("    }\n")
        f.write("}\n\n")

        # Section C: Entity Cleanup
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section C: Material & Property Database Purge\n")
        f.write("# ---------------------------------------------------------\n")
        f.write("*clearmark materials 1\ncatch { *createmark materials 1 all }\ncatch { *deletemark materials 1 }\n\n")
        f.write("*clearmark properties 1\ncatch { *createmark properties 1 all }\ncatch { *deletemark properties 1 }\n\n")

        # Section D: Create MAT8 Materials
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section D: Create Orthotropic Composite Materials (MAT8)\n")
        f.write("# ---------------------------------------------------------\n")
        for mat_id, mat in materials.items():
            f.write(f'*createsolverkeyword "MAT8" materials "Main Model"\n')
            f.write(f"*setvalue mats id={mat_id} include=0\n")
            f.write(f'*setvalue mats id={mat_id} name="{mat["name"]}"\n')
            f.write(f"*setvalue mats id={mat_id} STATUS=2 196={mat['E1']}\n")    # E1
            f.write(f"*setvalue mats id={mat_id} STATUS=2 197={mat['E2']}\n")    # E2
            f.write(f"*setvalue mats id={mat_id} STATUS=1 198={mat['NU12']}\n")  # NU12
            f.write(f"*setvalue mats id={mat_id} STATUS=2 199={mat['G12']}\n")   # G12
            f.write(f"*setvalue mats id={mat_id} STATUS=2 202={mat['RHO']}\n\n") # Density

        # Section E: Create PCOMP Property Cards and Assign to Components
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section E: Create & Bind Composite Property Cards (PCOMP)\n")
        f.write("# ---------------------------------------------------------\n")
        for idx, zone in enumerate(zones):
            prop_counter = idx + 1
            num_plies = len(zone["stack"])
            total_thickness = round(num_plies * zone["thickness"], 6)
            
            # Mid-plane offset calculation: Z0 = -Total_Thickness / 2
            z0_offset = round(-total_thickness / 2.0, 6)
            
            f.write(f'*createsolverkeyword "PCOMP" properties\n')
            f.write(f"*createmark props 1 -1\n")
            f.write(f"set prop_id_{prop_counter} [lindex [hm_getmark props 1] 0]\n")
            f.write(f'*setvalue props id=$prop_id_{prop_counter} name="{zone["name"]}"\n')
            f.write(f"*setvalue props id=$prop_id_{prop_counter} STATUS=2 3027={num_plies}\n")   # Ply Count
            f.write(f"*setvalue props id=$prop_id_{prop_counter} STATUS=2 3014={z0_offset}\n")   # Z0 Offset

            mat_array = " ".join([str(zone["material_id"])] * num_plies)
            thick_array = " ".join([str(zone["thickness"])] * num_plies)
            angle_array = " ".join([str(angle) for angle in zone["stack"]])
            
            f.write(f"*setvalue props id=$prop_id_{prop_counter} STATUS=2 3023={{mats {mat_array}}}\n")
            f.write(f"*setvalue props id=$prop_id_{prop_counter} STATUS=2 3024={{{thick_array}}}\n")
            f.write(f"*setvalue props id=$prop_id_{prop_counter} STATUS=2 3025={{{angle_array}}}\n")

            # Enable stress output (SOUT=YES) per ply
            for i in range(num_plies):
                f.write(f"*setvalue props id=$prop_id_{prop_counter} STATUS=2 ROW={i} 3026={{YES}}\n")

            # Bind property card to component
            f.write(f'*createmark components 1 "{zone["component"]}"\n')
            f.write(f"set comp_id_{prop_counter} [lindex [hm_getmark components 1] 0]\n")
            f.write(f"*setvalue comps id=$comp_id_{prop_counter} propertyid={{props $prop_id_{prop_counter}}}\n\n")

        # Section F: Element Card Image Formulations
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section F: Assign Shell Element Formulations (Card Images)\n")
        f.write("# ---------------------------------------------------------\n")
        for zone in zones:
            if "mesh_type" in zone:
                f.write(f'*createmark elems 1 "by comp name" "{zone["component"]}"\n')
                f.write(f'if {{[hm_marklength elems 1] > 0}} {{\n')
                f.write(f'    catch {{ *setvalue elems mark=1 cardimage="{zone["mesh_type"]}" }}\n')
                f.write(f'}}\n\n')

        # Section G: Boundary Conditions and Load Collectors
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section G: Setup Boundary Conditions (SPC) & Forces\n")
        f.write("# ---------------------------------------------------------\n")

        # Clean existing loadsteps and load collectors
        f.write("# Purge pre-existing load steps and load collectors\n")
        f.write("catch { *createmark loadsteps 1 \"Linear_Static\" }\n")
        f.write("catch { *deletemark loadsteps 1 }\n")
        f.write("catch { *createmark loadcols 1 \"SPC\" \"Forces\" }\n")
        f.write("catch { *deletemark loadcols 1 }\n\n")

        # Create Load Collector for SPC (Single Point Constraints)
        f.write("# Create SPC Load Collector\n")
        f.write("*createentity loadcols name=\"SPC\" color=4\n")
        for bc in bcs:
            f.write(f"# Apply Boundary Condition: {bc['name']}\n")
            f.write(f"*createmark nodes 1 \"by box\" {bc['x_min']} {bc['y_min']} {bc['z_min']} {bc['x_max']} {bc['y_max']} {bc['z_max']} 0 inside 0 1 0\n")
            f.write("if {[hm_marklength nodes 1] > 0} {\n")
            f.write(f"    *loadcreate 1 3 1 {bc['dofs_values']}\n")
            f.write(f"    puts \"[INFO] Boundary condition '{bc['name']}' applied successfully.\"\n")
            f.write("} else {\n")
            f.write(f"    puts \"[WARNING] No nodes found for BC '{bc['name']}'\"\n")
            f.write("}\n\n")

        # Create Load Collector for Forces
        f.write("# Create Forces Load Collector\n")
        f.write("*createentity loadcols name=\"Forces\" color=2\n")
        for load in loads:
            f.write(f"# Apply Point Load Vector: {load['name']}\n")
            f.write(f"*createmark nodes 1 \"by box\" {load['x_min']} {load['y_min']} {load['z_min']} {load['x_max']} {load['y_max']} {load['z_max']} 0 inside 0 1 0\n")
            f.write("if {[hm_marklength nodes 1] > 0} {\n")
            f.write(f"    *loadcreate 1 1 1 {load['vector']}\n")
            f.write(f"    puts \"[INFO] Force vector '{load['name']}' applied successfully.\"\n")
            f.write("} else {\n")
            f.write(f"    puts \"[WARNING] No nodes found for Force '{load['name']}'\"\n")
            f.write("}\n\n")

        # Section H: Create Linear Static Subcase (Load Step)
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section H: Create Linear Static Analysis Subcase\n")
        f.write("# ---------------------------------------------------------\n")
        f.write("*createentity loadsteps name=\"Linear_Static\"\n")
        f.write("*createmark loadsteps 1 \"Linear_Static\"\n")
        f.write("set step_id [lindex [hm_getmark loadsteps 1] 0]\n")
        f.write("*createmark loadcols 1 \"SPC\"\n")
        f.write("set spc_id [lindex [hm_getmark loadcols 1] 0]\n")
        f.write("*createmark loadcols 1 \"Forces\"\n")
        f.write("set force_id [lindex [hm_getmark loadcols 1] 0]\n")
        
        # Link SPC and Load collectors to the loadstep card attributes
        f.write("catch {*setvalue loadsteps id=$step_id STATUS=2 4149={loadcols $spc_id}}\n")
        f.write("catch {*setvalue loadsteps id=$step_id STATUS=2 4150={loadcols $force_id}}\n")
        f.write("catch {*setvalue loadsteps id=$step_id STATUS=2 4058=1}\n") # Linear Static analysis type
        f.write("puts \"[INFO] Load step 'Linear_Static' configured successfully.\"\n\n")

        # Section I: Geometric Verification & Visual Display Options
        f.write("# ---------------------------------------------------------\n")
        f.write("# Section I: Bounding Box Verification & Graphical Display\n")
        f.write("# ---------------------------------------------------------\n")
        for zone in zones:
            f.write(f'*createmark elems 1 "by comp name" "{zone["component"]}"\n')
            f.write(f'set num_elems [hm_marklength elems 1]\n')
            f.write(f'if {{$num_elems > 0}} {{\n')
            f.write(f'    set bbox [hm_getboundingbox elems 1 0 0 0]\n')
            f.write(f'    set x_min [lindex $bbox 0]\n')
            f.write(f'    set x_max [lindex $bbox 3]\n')
            f.write(f'    set y_min [lindex $bbox 1]\n')
            f.write(f'    set y_max [lindex $bbox 4]\n')
            f.write(f'    set x_dim [expr $x_max - $x_min]\n')
            f.write(f'    set y_dim [expr $y_max - $y_min]\n')
            f.write(f'    puts "[VERIFY] {zone["component"]}: ${{x_dim}}x${{y_dim}} mm ($num_elems elements)"\n')
            f.write(f'}} else {{\n')
            f.write(f'    puts "[WARNING] {zone["component"]}: NO ELEMENTS FOUND"\n')
            f.write(f'}}\n\n')

        # Enable 3D Composite Layer visualization representation
        f.write("*setoption show_composite_layers=2\n")
        f.write("catch {*createmark elems 1 all}\n")
        f.write("catch {*fitentitiestoscreen elems 1}\n\n")
        
        f.write("hm_answernext no\n")
        f.write('puts "\\n==================================================="\n')
        f.write('puts "[SUCCESS] MULTI-PATCH TCL SCRIPT COMPLETED"\n')
        f.write('puts "==================================================="\n')

    # -------------------------------------------------------------------------
    # 6. CONSOLE SUMMARY & HYPERMESH LAUNCH
    # -------------------------------------------------------------------------
    print(f"\n[SUCCESS] Script generated successfully: {OUTPUT_FILE}")
    print(f"  - Configured Zones   : {len(zones)} Quadrants")
    print(f"  - Materials Created  : {len(materials)} MAT8 cards")
    print(f"  - Analysis Subcase   : Linear_Static (SPC + Forces)")

    # Execute HyperMesh if path is valid
    if os.path.exists(HM_EXE):
        print("\n[INFO] Launching Altair HyperMesh...")
        subprocess.Popen([HM_EXE, "-tcl", OUTPUT_FILE])
    else:
        print(f"\n[ERROR] HyperMesh executable not found at: {HM_EXE}")

if __name__ == "__main__":
    generate_multi_patch_final()
