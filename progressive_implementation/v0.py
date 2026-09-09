import os
import subprocess

def generate_v0():
    """
    Generates the 'v0.tcl' automation script for Altair HyperMesh.
    
    Workflow steps:
    1. Define paths and structural laminate parameters.
    2. Write TCL commands to clear old data and import model file.
    3. Create MAT8 orthotropic material for Carbon Fiber Reinforced Polymer (CFRP).
    4. Create PCOMP composite property card and define 12-ply layup sequence.
    5. Assign composite property card to target structural component ('Zone_0').
    """

    # -------------------------------------------------------------------------
    # 1. PATH & PARAMETER CONFIGURATION
    # -------------------------------------------------------------------------
    # Path to HyperMesh executable (Altair 2025.1 Win64)
    HM_EXE = r"C:\Program Files\Altair\2025.1\hwdesktop\hm\bin\win64\hmopengl.exe"
    
    # Target HyperMesh database file (.hm)
    MODEL_PATH = r"c:\Users\Hugo\OneDrive\Escritorio\v0_1patch\v0.hm"
    
    # Output TCL script filename
    OUTPUT_FILE = "v0.tcl"
    
    # Composite Material & Property Parameters
    MAT_ID = 1                # Unique Material ID in OptiStruct
    THICKNESS = 0.15          # Individual ply thickness in mm
    PROP_NAME = "P_Single_Patch"  # Name for the PCOMP property card
    COMP_NAME = "Zone_0"      # Target component name in the model
    
    # Stacking Sequence: [25, -65, -25, 65] repeated 3 times -> 12 plies total
    stack = [25.0, -65.0, -25.0, 65.0] * 3  
    
    # Sanitize Windows file paths for TCL compatibility (forward slashes)
    tcl_model_path = MODEL_PATH.replace("\\", "/")

    print(f"[INFO] Generating automation TCL script: {OUTPUT_FILE}...")

    # -------------------------------------------------------------------------
    # 2. TCL SCRIPT FILE GENERATION
    # -------------------------------------------------------------------------
    with open(OUTPUT_FILE, "w") as f:
        # Script Header & Dialog Confirmation Suppression
        f.write("# HyperMesh TCL Automation Script - Generated via Python\n")
        f.write("hm_answernext yes\n\n")
        
        # Step 1: Read the base HyperMesh model database
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 1: Load Base HyperMesh Model\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f'*readfile "{tcl_model_path}" 0\n\n')

        # Step 2: Clean pre-existing entities to avoid naming/ID conflicts
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 2: Safe Cleanup of Existing Materials and Properties\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f"*clearmark materials 1\n")
        f.write(f'catch {{ *createmark materials 1 "Carbon_Auto" }}\n')
        f.write(f"catch {{ *deletemark materials 1 }}\n")
        f.write(f"*clearmark properties 1\n")
        f.write(f'catch {{ *createmark properties 1 "{PROP_NAME}" }}\n')
        f.write(f"catch {{ *deletemark properties 1 }}\n\n")

        # Step 3: Define MAT8 (2D Orthotropic Material for Composites)
        # OptiStruct Card Attributes:
        #   196: E1   (Longitudinal Modulus = 140,000 MPa)
        #   197: E2   (Transverse Modulus   = 10,000 MPa)
        #   198: NU12 (Poisson's Ratio      = 0.3)
        #   199: G12  (In-plane Shear Modulus = 5,000 MPa)
        #   191: RHO  (Material Density     = 1.6e-9 tonnes/mm^3)
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 3: Create MAT8 Orthotropic Composite Material\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f'*createsolverkeyword "MAT8" materials "Main Model"\n')
        f.write(f"*setvalue mats id={MAT_ID} include=0\n")
        f.write(f'*setvalue mats id={MAT_ID} name="Carbon_Auto"\n')
        f.write(f"*setvalue mats id={MAT_ID} STATUS=2 196=140000\n")  # E1
        f.write(f"*setvalue mats id={MAT_ID} STATUS=2 197=10000\n")   # E2
        f.write(f"*setvalue mats id={MAT_ID} STATUS=1 198=0.3\n")     # NU12
        f.write(f"*setvalue mats id={MAT_ID} STATUS=2 199=5000\n")    # G12
        f.write(f"*setvalue mats id={MAT_ID} STATUS=2 191=1.6e-9\n\n") # Density

        # Step 4: Create PCOMP Property Card
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 4: Create PCOMP Composite Shell Property Card\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f'*createsolverkeyword "PCOMP" properties\n')
        f.write(f"*createmark props 1 -1\n")  # Mark the last created property
        f.write(f"set prop_list [hm_getmark props 1]\n")
        f.write(f"set prop_id [lindex $prop_list 0]\n")
        f.write(f'*setvalue props id=$prop_id name="{PROP_NAME}"\n')
        f.write("puts \"[INFO] Created Property ID: $prop_id\"\n\n")

        # Step 5: Configure total ply count in PCOMP card (OptiStruct ID 3027)
        num_plies = len(stack)
        f.write("# ---------------------------------------------------------\n")
        f.write(f"# Step 5: Define Laminate Layup ({num_plies} Plies)\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f"*setvalue props id=$prop_id STATUS=2 3027={num_plies}\n\n")

        # Step 6: Assign arrays for Ply Materials, Thicknesses, and Angles
        # OptiStruct Card Attributes:
        #   3023: Material ID list per ply
        #   3024: Thickness list per ply
        #   3025: Orientation angle list per ply
        #   3026: Stress output flag (SOUT = YES)
        mat_array = " ".join([str(MAT_ID)] * num_plies)
        f.write(f"# Assign Material ID {MAT_ID} to all plies\n")
        f.write(f"*setvalue props id=$prop_id STATUS=2 3023={{mats {mat_array}}}\n\n")

        thick_array = " ".join([str(THICKNESS)] * num_plies)
        f.write(f"# Assign Ply Thickness ({THICKNESS} mm) to all plies\n")
        f.write(f"*setvalue props id=$prop_id STATUS=2 3024={{{thick_array}}}\n\n")

        angle_array = " ".join([str(angle) for angle in stack])
        f.write(f"# Assign Stacking Sequence Angles (deg)\n")
        f.write(f"*setvalue props id=$prop_id STATUS=2 3025={{{angle_array}}}\n\n")

        # Set Stress Output Flags (SOUT = YES for each ply row)
        f.write(f"# Enable Stress Output (SOUT=YES) across all plies\n")
        for i in range(num_plies):
            f.write(f"*setvalue props id=$prop_id STATUS=2 ROW={i} 3026={{YES}}\n")
        f.write("\n")

        # Step 7: Assign Composite Property Card to Target Component
        f.write("# ---------------------------------------------------------\n")
        f.write(f"# Step 7: Assign Property to Component ({COMP_NAME})\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f'*createmark components 1 "{COMP_NAME}"\n')
        f.write(f"set comp_list [hm_getmark components 1]\n")
        f.write(f"set comp_id [lindex $comp_list 0]\n")
        f.write("puts \"[INFO] Target Component ID: $comp_id\"\n\n")
        
        # Bind the PCOMP property card to the component ID
        f.write(f"*setvalue comps id=$comp_id propertyid={{props $prop_id}}\n\n")

        # Step 8: Verify Property Assignment in Log Output
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 8: Verify Property Assignment\n")
        f.write("# ---------------------------------------------------------\n")
        f.write('set assigned_prop [hm_getvalue comps id=$comp_id dataname=propertyid]\n')
        f.write('puts "[VERIFY] Property ID assigned to component: $assigned_prop"\n\n')

        # Step 9: Visual Highlight and Script Completion Wrap-up
        f.write("# ---------------------------------------------------------\n")
        f.write("# Step 9: Highlight Component & Finalize\n")
        f.write("# ---------------------------------------------------------\n")
        f.write(f'*createmark components 1 "{COMP_NAME}"\n')
        f.write("hm_highlightmark components 1 h\n")
        f.write("hm_answernext no\n")
        f.write("puts \"===================================================\"\n")
        f.write("puts \"[SUCCESS] TCL AUTOMATION SCRIPT EXECUTED SUCCESSFULLY\"\n")
        f.write("puts \"===================================================\"\n")

    # -------------------------------------------------------------------------
    # 3. CONSOLE SUMMARY & HYPERMESH LAUNCH
    # -------------------------------------------------------------------------
    print(f"[SUCCESS] TCL script successfully generated: {OUTPUT_FILE}")
    print(f"  - Material ID       : {MAT_ID} (Carbon_Auto)")
    print(f"  - Property Name     : {PROP_NAME}")
    print(f"  - Component Name    : {COMP_NAME}")
    print(f"  - Total Plies       : {num_plies}")
    print(f"  - Stacking Sequence : {stack}")

    # Launch HyperMesh in non-blocking background process
    if os.path.exists(HM_EXE):
        print("\n[INFO] Launching Altair HyperMesh...")
        subprocess.Popen([HM_EXE, "-tcl", OUTPUT_FILE])
    else:
        print(f"\n[ERROR] Executable path not found: {HM_EXE}")

if __name__ == "__main__":
    generate_v0()
