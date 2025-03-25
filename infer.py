import os
from fastapi import FastAPI, File, UploadFile, Form, HTTPException
from fastapi.responses import FileResponse

import uvicorn
from rdkit import Chem
import MDAnalysis as mda
from inference import main as inference_main, get_parser as inference_get_parser

app = FastAPI()

def run_inference(pdb_path: str, ligand: str, complex_name: str):
    output_dir = "./results"
    os.makedirs(output_dir, exist_ok=True)

    parser = inference_get_parser()
    args = parser.parse_args([])
    
    args.protein_path = pdb_path
    args.ligand_description = ligand
    args.complex_name = complex_name
    args.save_visualisation = False
    args.out_dir = output_dir

    inference_main(args)

    return os.path.join(output_dir, complex_name, "rank1.sdf")

@app.post("/merge")
async def merge(
    pdb_file: UploadFile = File(...),
    ligand_file: UploadFile = File(...),
    output_path: str = Form(...)
):
    """
    Merge a protein and a ligand PDB file and save the result to a new file.

    Args:
        pdb_file: The protein PDB file.
        ligand_file: The ligand PDB file.
        output_path: The path to save the merged file.
    """
    try:
        ligand = Chem.MolFromMolFile(ligand_file.file)
        Chem.MolToPDBFile(ligand, output_path)

        protein = mda.Universe(pdb_file.file)
        ligand = mda.Universe(ligand_file.file)
        combined = mda.Merge(protein.atoms, ligand.atoms)
        combined.atoms.write(output_path)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/infer")
async def infer(
    pdb_file: UploadFile = File(...),  # Accepts a file
    ligand: str = Form(...),  # Accepts form data
    complex_name: str = Form(...),
    mock: bool = Form(False)
):
    """
    Run a DiffDock inference on a protein and a ligand.
    Args:
        pdb_file (UploadFile, optional): The protein PDB file. Defaults to File(...).
        ligand (str, optional): The ligand PDB file. Defaults to Form(...).
        complex_name (str, optional): The name of the complex. Defaults to Form(...).
        mock (bool, optional): Whether to run a mock inference. Defaults to Form(False).

    Raises:
        HTTPException: _description_
        HTTPException: _description_

    Returns:
        FileResponse: The predicted complex.
    """
    if mock:
        return "mock_results"

    input_dir = "test_folding"
    os.makedirs(input_dir, exist_ok=True)

    pdb_path = os.path.join(input_dir, f"{complex_name}.pdb")

    # Save uploaded file
    with open(pdb_path, "wb") as f:
        f.write(await pdb_file.read())

    try:
        output_path = run_inference(pdb_path, ligand, complex_name)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
    if not os.path.exists(output_path):
        raise HTTPException(status_code=404, detail="Output file not found.")

    return FileResponse(output_path, media_type="chemical/x-sdf", filename=f"{complex_name}_result.sdf")

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)