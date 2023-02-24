
#!/usr/bin/env python
"""
An animated image
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import types
fig = plt.figure()
def f(x, y):
    return np.sin(x) + np.cos(y)
x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
ims = []
for i in range(60):
    x += np.pi / 15.
    y += np.pi / 20.
    im = plt.contour(f(x, y))
    def setvisible(self,vis):
        print(self)
        for c in self[0].collections: c.set_visible(self[1])
    im.set_visible = types.MethodType(setvisible,(im,None,))
    # im.set_visible = setvisible(im,None)
    im.axes = plt.gca()
    ims.append([im])
# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
#     repeat_delay=1000)
# ani.save('dynamic_images.mp4')

ani = animation.FuncAnimation(fig,ims ,frames = 60, blit=False,repeat=False)
writer = animation.writers['ffmpeg'](fps=6)
ani.save("zer.py",writer=writer,dpi=200)

plt.show()
