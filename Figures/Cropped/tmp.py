# %%

from PIL import Image
import os

def crop_image(input_image_path, output_dir, crop_areas):
    """
    Crop the input image into parts defined by crop_areas and save them.

    Parameters:
    - input_image_path: Path to the input image.
    - output_dir: Directory where cropped images will be saved.
    - crop_areas: Dictionary defining crop areas for each part.
                  Each key is the name of the part (e.g., 'A', 'B') and each value is a tuple (left, upper, right, lower).
    """
    image = Image.open(input_image_path)
    
    os.makedirs(output_dir, exist_ok=True)
    
    for part, area in crop_areas.items():
        cropped_image = image.crop(area)
        cropped_image.save(os.path.join(output_dir, f"Supplementary_Figure_8{part}.png"), format='PNG')
        cropped_image.save(os.path.join(output_dir, f"Supplementary_Figure_8{part}.tiff"), format='TIFF')

# Define the input image path and output directory
input_image_path = "/BiO/Access/kyungwhan1998/transcriptomic_clock/Figures/Supplementary_Figure_8.png"
output_dir = "/BiO/Access/kyungwhan1998/transcriptomic_clock/Figures/Cropped"

image = Image.open(input_image_path)

## Supplementary_Figure_8
crop_areas = {
    'A': (0, 0, image.size[0]/2+250, image.size[1]/3-50),   # Coordinates for part A
    'B': (0, image.size[1]/3-50, image.size[0]/2+250, 2*image.size[1]/3-100),  # Coordinates for part B
    'C': (0, 2*image.size[1]/3-100, image.size[0]/2+250, image.size[1]),  # Coordinates for part C
    'D': (image.size[0]/2+250, 0, image.size[0], image.size[1]),  # Coordinates for part D
}


# Crop the image and save the parts
crop_image(input_image_path, output_dir, crop_areas)
# %%
