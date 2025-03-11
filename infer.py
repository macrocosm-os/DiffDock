import os

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel

from inference import main as inference_main, get_parser as inference_get_parser

app = FastAPI()

class InferRequest(BaseModel):
    pdb_input_path: str
    ligand: str
    complex_name: str
    mock: bool

def run_inference(pdb_input_path: str, ligand: str, complex_name: str):
    # Ensure the input directory exists
    input_dir = "test_folding"
    os.makedirs(input_dir, exist_ok=True)

    # Ensure the output directory exists
    output_dir = "./results"
    os.makedirs(output_dir, exist_ok=True)
    
    # Get default arguments and modify them
    parser = inference_get_parser()
    args = parser.parse_args([])  # Start with empty args
    
    # Set up the necessary arguments
    args.protein_path = pdb_input_path
    args.ligand_description = ligand
    args.complex_name = complex_name
    args.save_visualisation = False
    args.out_dir = output_dir
    
    # Call the inference main function directly
    inference_main(args)
    
    # Return the path to the output file
    return os.path.join(output_dir, complex_name, "rank1.sdf")

def run_inference_mock():
    return "mock_results"

@app.post("/infer")
async def infer(request: InferRequest):
    # return {"message": "Request received!"}

    if request.mock:
        return run_inference_mock()

    try:
        output_path = run_inference(request.pdb_input_path, request.ligand, request.complex_name)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
    if not os.path.exists(output_path):
        raise HTTPException(status_code=404, detail="Output file not found.")
    
    # Return the file as a response with the appropriate MIME type
    return FileResponse(output_path, media_type="chemical/x-sdf", filename=f"{request.complex_name}_result.sdf")
