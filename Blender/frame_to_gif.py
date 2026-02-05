import sys
import argparse
import os
import re
from pathlib import Path
from PIL import Image
import cv2

# Get user input for output type
output_type = input("Enter the type of the output (gif, mp4): ")
starting_frame = 0

compression = round(float(input("compression (0-1): ")),1 )
print(compression)
rad = int(input("Radius: "))
print(rad)

# Parse command-line arguments
parser = argparse.ArgumentParser()
argv = sys.argv
if len(argv) > 1:
    filename = argv[1]
else:
    raise ValueError("Please provide a filename.")

# Change directory to images folder
os.chdir(filename[:-4] + '/images')

# Natural sorting function
def natural_sort_key(filename):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', filename)]

# Get PNG files and sort them
files = sorted([f for f in os.listdir("./") if f.endswith(".png")], key=natural_sort_key)

if output_type == 'movie':
    first_image = cv2.imread(files[0])
    height, width, _ = first_image.shape

    upscale_factor = 2
    new_width, new_height = width * upscale_factor, height * upscale_factor

    output_video = f"./animated_{filename[:-4]}.mp4"
    if os.path.exists(output_video):
        os.remove(output_video)

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    fps = 10
    video_writer = cv2.VideoWriter(output_video, fourcc, fps, (new_width, new_height))

    for idx, file in enumerate(files[starting_frame:], start=starting_frame):
        frame = cv2.imread(file)
        frame = cv2.resize(frame, (new_width, new_height), interpolation=cv2.INTER_CUBIC)
        
        # Add text in the top-left corner
        text1 = f"c = {compression}"
        text2 = f"R = {rad} micro meters"
        text3 = f"Frame: {idx + 1}"
        font = cv2.FONT_HERSHEY_TRIPLEX  # Closest to Times Roman in OpenCV
        font_scale = 2
        font_color = (0, 0, 0)  # Black
        thickness = 2
        margin = 20
        
        text_x1 = margin
        text_y1 = margin + cv2.getTextSize(text1, font, font_scale, thickness)[0][1]
        text_x2 = margin
        text_y2 = text_y1 + cv2.getTextSize(text2, font, font_scale, thickness)[0][1] + margin
        text_x3 = margin
        text_y3 = text_y2 + cv2.getTextSize(text3, font, font_scale, thickness)[0][1] + margin

        cv2.putText(frame, text1, (text_x1, text_y1), font, font_scale, font_color, thickness)
        cv2.putText(frame, text2, (text_x2, text_y2), font, font_scale, font_color, thickness)
        cv2.putText(frame, text3, (text_x3, text_y3), font, font_scale, font_color, thickness)

        video_writer.write(frame)
    video_writer.release()
    print(f"Video saved as {output_video}")

else:  # GIF creation
    images = [Image.open(file) for file in files[starting_frame:]]

    # Reduce colors to optimize GIF size
    images = [img.convert("P", palette=Image.ADAPTIVE, colors=128) for img in images]

    output_gif = f"./animated_{filename[:-4]}.gif"
    if os.path.exists(output_gif):
        os.remove(output_gif)

    # Save GIF with optimization
    images[0].save(
        output_gif,
        save_all=True,
        append_images=images[1:],
        optimize=True,
        duration=150,  # Slightly slower to reduce size
        loop=0
    )

    print(f"GIF saved as {output_gif}, final size: {os.path.getsize(output_gif) / (1024*1024):.2f} MB")
