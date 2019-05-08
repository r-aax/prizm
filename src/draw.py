# -*- coding: utf-8 -*-
"""
Drawing examples.

Created on Tue May  7 09:58:45 2019

@author: Rybakov
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import aggdraw
from PIL import Image, ImageDraw, ImageFont

#---------------------------------------------------------------------------------------------------
# Other functions.
#---------------------------------------------------------------------------------------------------

def test_pil():
    """
    PIL test.
    """

    # Create image and drawer.
    img = Image.new('RGB', (600, 300), color = (73, 109, 137))
    c = ImageDraw.Draw(img)

    # Draw text.
    fnt = ImageFont.truetype('C:\Windows\Fonts\Arial.ttf', 15)
    c.text((10, 10), "Arial", font = fnt, fill = (255, 255, 0))
    #
    fnt = ImageFont.truetype('C:\Windows\Fonts\lucon.ttf', 18)
    c.text((200, 80), "Lucida Console", font = fnt, fill = (0, 255, 255))
    #
    fnt = ImageFont.truetype('C:\Windows\Fonts\COLONNA.ttf', 20)
    c.text((400, 170), "Colonna MT", font = fnt, fill = (255, 0, 0))
    #
    fnt = ImageFont.truetype('C:\Windows\Fonts\OLDENGL.ttf', 14)
    c.text((350, 250), "Old English Text", font = fnt, fill = (60, 70, 80))

    # Draw primitives.
    c.ellipse((50, 50, 200, 100), fill = 'red', outline = 'blue')
    c.line((10, 10, 590, 290))

    # Save and show.
    img.save('test.png')
    img.show()

#---------------------------------------------------------------------------------------------------

def test_aggdraw():
    """
    aggdraw test.
    """

    # Create image and drawer.
    img = Image.new('RGB', (600, 300), color = (73, 109, 137))
    c = aggdraw.Draw(img)

    c.setantialias(True)

    # Draw aggdraw primitives.
    c.arc((500, 200, 550, 280), 85, 320, aggdraw.Pen('red', 2.1))
    c.chord((400, 200, 450, 290), 10, 290, aggdraw.Pen('orange', 2.5), aggdraw.Brush('blue'))
    c.ellipse((50, 50, 200, 100), aggdraw.Pen('orange', 3.5), aggdraw.Brush('green'))
    c.line((10, 10, 590, 290), aggdraw.Pen('blue', 2.9))
    c.line((5, 5, 100, 20, 30, 40, 400, 250, 550, 250, 300, 100), aggdraw.Pen('red', 0.5))
    #
    path = aggdraw.Path()
    path.moveto(0, 0)
    path.curveto(200, 250, 500, 50, 400, 250)
    path.curveto(0, 60, 40, 100, 100, 100)
    c.path(path, aggdraw.Pen('pink', 5.0))
    #
    c.pieslice((100, 100, 300, 250), 50, 350,
               aggdraw.Pen('green', 2.5), aggdraw.Brush('yellow'))
    c.polygon((590, 10, 590, 290, 580, 290, 580, 20, 500, 10),
              aggdraw.Pen('blue', 1.0), aggdraw.Brush('pink'))
    c.rectangle((400, 200, 500, 300), aggdraw.Pen('black', 2.8), aggdraw.Brush('orange'))

    # Flush
    c.flush()

    # Save and show.
    img.save('test.png')
    img.show()

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    #test_pil()
    test_aggdraw()

#---------------------------------------------------------------------------------------------------