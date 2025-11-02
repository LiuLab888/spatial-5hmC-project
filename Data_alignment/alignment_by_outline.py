
#========= Load libraries ========
# Import necessary libraries
import numpy as np
from skimage import io, color
from skimage.measure import block_reduce
from skimage.morphology import binary_closing, disk
import SimpleITK as sitk
import matplotlib.pyplot as plt
import pandas as pd

#========= Load images ========
# 1. Read in the images and convert to grayscale
img_meth = io.imread("MethDotplot.jpg")  # Hydroxymethylation image
img_rna  = io.imread("RnaDotplot.jpg")   # RNA image

# Convert images to grayscale
gray_meth = color.rgb2gray(img_meth)
gray_rna = color.rgb2gray(img_rna)

# Display original grayscale images side by side
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].imshow(gray_meth, cmap="gray")  # Display Hydroxymethylation image
ax[0].set_title("Original Hydroxymethylation (Gray)")
ax[1].imshow(gray_rna, cmap="gray")   # Display RNA image
ax[1].set_title("Original RNA (Gray)")
plt.show()

# 2. Downsize images to 192x192 pixels
# Calculate block size for downsampling (192x192 target resolution)
block_size_y = gray_meth.shape[0] // 192
block_size_x = gray_meth.shape[1] // 192

# Downsample both images using block_reduce (average block value)
meth_192 = block_reduce(gray_meth, block_size=(block_size_y, block_size_x), func=np.mean)
rna_192 = block_reduce(gray_rna, block_size=(block_size_y, block_size_x), func=np.mean)

# 2.1 Visualization of downsampled images
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].imshow(meth_192, cmap="gray")  # Display downsampled Hydroxymethylation image
ax[0].set_title("Downsampled 192x192 Hydroxymethylation")
ax[1].imshow(rna_192, cmap="gray")   # Display downsampled RNA image
ax[1].set_title("Downsampled 192x192 RNA")
plt.show()

#========= Emphasize outline ========
# 3. Perform binary edge closing operation
edges_meth = binary_closing(meth_192 > 0.999, disk(1))  # Apply closing operation to Hydroxymethylation image
edges_rna = binary_closing(rna_192 > 0.999, disk(1))    # Apply closing operation to RNA image

# 3.1 Visualization of the edge-processed images
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].imshow(edges_meth, cmap="gray")  # Display binary edge of Hydroxymethylation image
ax[0].set_title("Meth 192x192 Edge (binary)")
ax[1].imshow(edges_rna, cmap="gray")   # Display binary edge of RNA image
ax[1].set_title("RNA 192x192 Edge (binary)")
plt.show()

# 4. Convert images to SimpleITK images for registration
fixed_img = sitk.GetImageFromArray(edges_rna.astype(np.uint8))  # Convert RNA edge image to SimpleITK format
moving_img = sitk.GetImageFromArray(edges_meth.astype(np.uint8))  # Convert Hydroxymethylation edge image to SimpleITK format

# 5. Perform image registration using affine transformation
elastix = sitk.ElastixImageFilter()  # Create elastix image filter for registration
elastix.SetFixedImage(fixed_img)  # Set fixed image (RNA)
elastix.SetMovingImage(moving_img)  # Set moving image (Hydroxymethylation)

# Set affine transformation parameters
param_map = sitk.GetDefaultParameterMap("affine")
param_map["Size"] = ['192', '192']
param_map["MaximumNumberOfIterations"] = ["1000"]
param_map["NumberOfResolutions"] = ["2"] #Test various values until optimal alignment
param_map["MaximumStepLength"] = ["2"] #Test various values until optimal alignment
param_map["AutomaticScalesEstimation"] = ["true"]
param_map["AutomaticTransformInitialization"] = ["true"]
elastix.SetParameterMap(param_map)
elastix.Execute()  # Run elastix to perform the registration

# Get transformation parameters after registration
transform_param_map = elastix.GetTransformParameterMap()

# 6. Apply the transformation to the original meth image using Transformix
moving_for_transform = sitk.GetImageFromArray((meth_192 * 255).astype(np.uint8))  # Prepare image for transformation

# Apply the transformation to the moving image (Hydroxymethylation)
transformix = sitk.TransformixImageFilter()
transformix.SetMovingImage(moving_for_transform)
transformix.SetTransformParameterMap(transform_param_map)
transformix.Execute()  # Execute transformation

# Get the result image from the transformation and normalize it
registered = transformix.GetResultImage()
registered_np = sitk.GetArrayFromImage(registered) / 255.0  # Normalize the image (0-1)

# 7. Visualize the images before and after alignment
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
# Before registration: Overlay original RNA and Hydroxymethylation images
ax[0].imshow(rna_192, cmap="Reds", alpha=0.5)  # RNA image with red colormap
ax[0].imshow(meth_192, cmap="Blues", alpha=0.5)  # Hydroxymethylation image with blue colormap
ax[0].set_title("Before Alignment")

# After registration: Overlay registered Hydroxymethylation and RNA images
ax[1].imshow(rna_192, cmap="Reds", alpha=0.5)  # RNA image with red colormap
ax[1].imshow(registered_np, cmap="Blues", alpha=0.5)  # Registered Hydroxymethylation image with blue colormap
ax[1].set_title("After Alignment")
plt.tight_layout()
plt.show()

# 8. Extract coordinates for comparison (pre- and post-alignment)
# Generate coordinate grid for registered image
H, W = registered_np.shape
xx, yy = np.meshgrid(np.arange(W), np.arange(H))  # Create a grid of pixel coordinates
coords_grid = np.stack([xx.flatten(), yy.flatten()], axis=1)  # Flatten grid into coordinates (H*W, 2)

# Extract non-zero pixel positions in the registered image (transformed coordinates)
mask = registered_np.flatten() > 0
coords_trans = coords_grid[mask]  # Get transformed coordinates (non-zero pixels)

# Extract corresponding coordinates in the original RNA image (original coordinates)
rna_mask = rna_192.flatten() > 0
coords_orig = coords_grid[rna_mask]  # Get original coordinates (non-zero pixels in RNA)

# Ensure both coordinate sets have the same length for comparison (conservative approach)
n_points = min(len(coords_orig), len(coords_trans))  # Match the number of points
coords_orig = coords_orig[:n_points]
coords_trans = coords_trans[:n_points]

# 9. Save the coordinates to a CSV file for further analysis
df = pd.DataFrame({
    "X_orig": coords_orig[:, 0],  # Original X coordinates
    "Y_orig": coords_orig[:, 1],  # Original Y coordinates
    "X_trans": coords_trans[:, 0],  # Transformed X coordinates
    "Y_trans": coords_trans[:, 1]   # Transformed Y coordinates
})
df.to_csv("dot_movement_from_registered_np.csv", index=False)  # Save as CSV
