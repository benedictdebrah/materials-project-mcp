from dataclasses import dataclass




@dataclass
class ElectronicBandStructurePrompt:
    """
    Dataclass to store the prompt template for electronic band structure tool usage.
    Contains the comprehensive template with instructions, parameters, and examples.
    """
    template: str = """
        # Electronic Band Structure Tool Usage Guide

        ## Tool Overview
        The get_electronic_bandstructure tool generates electronic band structure plots for materials from the Materials Project database. It returns the plot as a base64-encoded PNG image within a JSON response structure.

        ## Tool Parameters
        - material_id (required): Materials Project ID (e.g., 'mp-149', 'mp-22526')
        - path_type (optional): K-point path type:
        - 'setyawan_curtarolo' (default) - Standard path for cubic systems
        - 'hinuma' - Standard path for hexagonal systems
        - 'latimer_munro' - Alternative path for cubic systems
        - 'uniform' - Uniform k-point sampling (not recommended for plotting)

        ## Expected Output Format
        The tool returns a JSON object with the following structure:
        
        {
        "success": true,
        "material_id": "mp-149",
        "image_base64": "iVBORw0KGgoAAAANSUhEUgAAB...[~500KB-2MB base64 string]",
        "metadata": {
            "path_type": "setyawan_curtarolo",
            "description": "Band structure plot for material mp-149 using setyawan_curtarolo path",
            "width": 1200,
            "height": 800
        }
        }
        

        Note: The image_base64 field contains a very long base64 string (typically 500KB-2MB). For brevity, examples show truncated versions.

        ## How to Parse and Display the Response

        ### Method 1: Extract and Display Base64 Image (Python)
        
        import json
        import base64
        from PIL import Image
        import io
        import matplotlib.pyplot as plt

        def display_bandstructure(response):
            # Parse the JSON response
            data = json.loads(response) if isinstance(response, str) else response
            
            if not data.get("success"):
                print("Error: Tool execution failed")
                return
            
            # Extract base64 image data (this will be a very long string)
            image_base64 = data["image_base64"]
            material_id = data["material_id"]
            metadata = data["metadata"]
            
            # Decode base64 to image
            image_bytes = base64.b64decode(image_base64)
            image = Image.open(io.BytesIO(image_bytes))
            
            # Display the image
            plt.figure(figsize=(12, 8))
            plt.imshow(image)
            plt.axis('off')
            plt.title(f"Band Structure: {material_id} ({metadata['path_type']})")
            plt.tight_layout()
            plt.show()
            
            print(f"Material ID: {material_id}")
            print(f"Path Type: {metadata['path_type']}")
            print(f"Image Size: {metadata['width']} × {metadata['height']}")
        
        """