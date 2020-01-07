import imageio

with imageio.get_writer('motion.gif', mode='I', duration=0.1) as writer:
    #for pi in parr:
    for pi in range(50):
        image = imageio.imread(f'shade_motion_{pi:03d}.png')
        writer.append_data(image)
