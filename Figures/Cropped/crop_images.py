# %%
import os

from PIL import Image


def crop_image(input_image_path, figure, output_dir, crop_areas):
    """
    Crop the input image into parts defined by crop_areas and save them.

    Parameters:
    - input_image_path: Path to the input image.
    - output_dir: Directory where cropped images will be saved.
    - crop_areas: Dictionary defining crop areas for each part.
    """
    image = Image.open(input_image_path)
    os.makedirs(output_dir, exist_ok=True)
    
    for part, area in crop_areas.items():
        cropped_image = image.crop(area)
        cropped_image.save(os.path.join(output_dir, f"{figure}{part}.png"), format='PNG')
        cropped_image.save(os.path.join(output_dir, f"{figure}{part}.tiff"), format='TIFF')

# %%
# Define the input image path and output directory
figure = "Figure_1"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Figure_1
crop_areas = {
    'A': (0, 0, image.size[0]/3-70, image.size[1]/3-20),   # Coordinates for part A
    'B': (0, image.size[1]/3, image.size[0]/3-70, 2*image.size[1]/3-40),  # Coordinates for part B
    'C': (0, 2*image.size[1]/3-40, image.size[0]/3-70, image.size[1]-60),  # Coordinates for part C
    'D': (image.size[0]/3-70, 0, 2*image.size[0]/3-280, image.size[1]),  # Coordinates for part D
    'E': (2*image.size[0]/3-280, 0, image.size[0], image.size[1]/2-120),  # Coordinates for part E
    'F': (2*image.size[0]/3-280, image.size[1]/2-120, image.size[0], image.size[1]),  # Coordinates for part F
}

crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Figure_2"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Figure_2
crop_areas = {
    'A': (0, 0, 2*image.size[0]/3-50, image.size[1]),   # Coordinates for part A
    'B': (2*image.size[0]/3-50, 0, image.size[0], image.size[1]),  # Coordinates for part B
}

crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Figure_3"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Figure_2
crop_areas = {
    'A': (0, 0, 2*image.size[0]/3-50, image.size[1]/2+30),   # Coordinates for part A
    'B': (2*image.size[0]/3-50, 0, image.size[0], image.size[1]/2+30),   # Coordinates for part B
    'C': (0, image.size[1]/2+30, 2*image.size[0]/3-50, image.size[1]),  # Coordinates for part C
    'D': (2*image.size[0]/3-50, image.size[1]/2+30, image.size[0], image.size[1]),  # Coordinates for part D
}
# each value is a tuple (left, upper, right, lower)
crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Supplementary_Figure_2"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_2
crop_areas = {
    'A': (0, 0, image.size[1]/2+20, image.size[0]/3-90),   # Coordinates for part A
    'B': (image.size[1]/2+20, 0, image.size[1]+20, image.size[0]/3-90),  # Coordinates for part B
    'C': (image.size[1]+20, 0, image.size[0], image.size[0]/3-90),  # Coordinates for part C
    'D': (0, image.size[0]/3-90, image.size[0], image.size[0]),  # Coordinates for part D
}
# each value is a tuple (left, upper, right, lower)
crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Supplementary_Figure_3"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_3
crop_areas = {
    'A': (0, 0, image.size[0]/3, image.size[1]/3),   # Coordinates for part A
    'B': (image.size[0]/3, 0, 2*image.size[0]/3, image.size[1]/3),  # Coordinates for part B
    'C': (2*image.size[0]/3, 0, image.size[0], image.size[1]/3),  # Coordinates for part C
    'D': (0, image.size[1]/3, image.size[0]/3, 2*image.size[1]/3),  # Coordinates for part D
    'E': (image.size[0]/3, image.size[1]/3, 2*image.size[0]/3, 2*image.size[1]/3),  # Coordinates for part E
    'F': (2*image.size[0]/3, image.size[1]/3, image.size[0], 2*image.size[1]/3),  # Coordinates for part F
    'G': (0, 2*image.size[1]/3, image.size[0]/3, image.size[1]),  # Coordinates for part G
    'H': (image.size[0]/3, 2*image.size[1]/3, 2*image.size[0]/3, image.size[1]),  # Coordinates for part H
    'I': (2*image.size[0]/3, 2*image.size[1]/3, image.size[0], image.size[1])  # Coordinates for part I
}
# each value is a tuple (left, upper, right, lower)
crop_image(input_image_path, figure, output_dir, crop_areas)
#%%
# Define the input image path and output directory
figure = "Supplementary_Figure_4"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_4
crop_areas = {
    'A': (0, 0, image.size[0]/2, image.size[1]),   # Coordinates for part A
    'B': (image.size[0]/2, 0, image.size[0], image.size[1]),  # Coordinates for part B
}
# each value is a tuple (left, upper, right, lower)
crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Supplementary_Figure_5"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_5
crop_areas = {
    'A': (0, 0, image.size[0]/2, image.size[1]/4+13),   # Coordinates for part A
    'B': (image.size[0]/2, 0, image.size[0], image.size[1]/4+13),  # Coordinates for part B
    'C': (0, image.size[1]/4+13, image.size[0]/2, image.size[1]/2+10),  # Coordinates for part C
    'D': (image.size[0]/2, image.size[1]/4+13, image.size[0], image.size[1]/2+5),  # Coordinates for part D
    'E': (0, image.size[1]/2+5, image.size[0]/2, 3*image.size[1]/4-15),  # Coordinates for part E
    'F': (image.size[0]/2, image.size[1]/2+5, image.size[0], 3*image.size[1]/4-15),  # Coordinates for part F
    'G': (0, 3*image.size[1]/4-15, image.size[0]/2, image.size[1]),  # Coordinates for part G
    'H': (image.size[0]/2, 3*image.size[1]/4-15, image.size[0], image.size[1]),  # Coordinates for part H
}

crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Supplementary_Figure_6"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_6
crop_areas = {
    'A': (0, 0, image.size[0]/2, image.size[1]/4),   # Coordinates for part A
    'B': (image.size[0]/2, 0, image.size[0], image.size[1]/4),  # Coordinates for part B
    'C': (0, image.size[1]/4, image.size[0]/2, image.size[1]/2),  # Coordinates for part C
    'D': (image.size[0]/2, image.size[1]/4, image.size[0], image.size[1]/2),  # Coordinates for part D
    'E': (0, image.size[1]/2, image.size[0]/2, 3*image.size[1]/4),  # Coordinates for part E
    'F': (image.size[0]/2, image.size[1]/2, image.size[0], 3*image.size[1]/4),  # Coordinates for part F
    'G': (0, 3*image.size[1]/4, image.size[0]/2, image.size[1])  # Coordinates for part G
}

crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Supplementary_Figure_8"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_8
crop_areas = {
    'A': (0, 0, image.size[0]/2+150, image.size[1]/3),   # Coordinates for part A
    'B': (0, image.size[1]/3, image.size[0]/2+150, 2*image.size[1]/3-30),  # Coordinates for part B
    'C': (0, 2*image.size[1]/3-30, image.size[0]/2+150, image.size[1]),  # Coordinates for part C
    'D': (image.size[0]/2+150, 0, image.size[0], image.size[1]),  # Coordinates for part D
}
crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Supplementary_Figure_10"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_8
crop_areas = {
    'A': (0, 0, image.size[0], image.size[1]/2),   # Coordinates for part A
    'B': (0, image.size[1]/2, image.size[0], image.size[1]),  # Coordinates for part B
}
crop_image(input_image_path, figure, output_dir, crop_areas)

# %%
# Define the input image path and output directory
figure = "Supplementary_Figure_11"
input_image_path = f"/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/{figure}.png"
output_dir = "/BiO/Access/kyungwhan1998/Backup/transcriptomic_clock/Figures/Cropped"
image = Image.open(input_image_path)
## Supplementary_Figure_8
crop_areas = {
    'A': (0, 0, image.size[0]/2, image.size[1]),   # Coordinates for part A
    'B': (image.size[0]/2, 0, image.size[0], image.size[1]),  # Coordinates for part B
}
crop_image(input_image_path, figure, output_dir, crop_areas)

# each value is a tuple (left, upper, right, lower)