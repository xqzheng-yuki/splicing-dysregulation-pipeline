# pip install pillow
import argparse
from PIL import Image
import os
import math

def get_arg_or_input(arg_value, prompt):
    return arg_value if arg_value else input(prompt)

parser = argparse.ArgumentParser()
parser.add_argument("--target", help="Your target")
args = parser.parse_args()

target = get_arg_or_input(args.target, "Enter your target file: ")

# target = "TPM_correlation"
image_folder = "/mnt/gtklab01/xiaoqing/result_for_thesis"
output_path = os.path.join(image_folder, f"{target}_combined.png")

# Get all relevant image files
image_files = sorted([
    os.path.join(image_folder, f) for f in os.listdir(image_folder)
    if f.startswith(target) and f.endswith(".png")
])

# Load images
images = [Image.open(f) for f in image_files]

# Step 1: Determine the max width and height
max_width = max(img.width for img in images)
max_height = max(img.height for img in images)

# Step 2: Scale each image proportionally to fit max dimensions
resized_images = []
for img in images:
    # Calculate scale factor to match either width or height
    scale_w = max_width / img.width
    scale_h = max_height / img.height
    scale = min(scale_w, scale_h)  # scale up proportionally without exceeding max

    new_size = (int(img.width * scale), int(img.height * scale))
    resized_img = img.resize(new_size, Image.ANTIALIAS)

    # Step 3: Center-pad the resized image to final dimensions
    final_img = Image.new("RGB", (max_width, max_height), "white")
    paste_x = (max_width - new_size[0]) // 2
    paste_y = (max_height - new_size[1]) // 2
    final_img.paste(resized_img, (paste_x, paste_y))
    resized_images.append(final_img)

# Step 4: Stitch into grid
cols = 4
rows = math.ceil(len(resized_images) / cols)

combined_image = Image.new("RGB", (cols * max_width, rows * max_height), "white")

for idx, img in enumerate(resized_images):
    row = idx // cols
    col = idx % cols
    x = col * max_width
    y = row * max_height
    combined_image.paste(img, (x, y))

# Save final combined image
combined_image.save(output_path)
print(f"✅ Combined image saved at: {output_path}")
