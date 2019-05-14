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
# Class drawer.
#---------------------------------------------------------------------------------------------------

class Drawer:
    """
    Drawer based on aggdraw.
    """

#---------------------------------------------------------------------------------------------------

    def __init__(self,
                 draw_area = (0.0, 0.0, 100.0, 100.0),
                 pic_size = (640, 480),
                 margins = (10, 10),
                 invert = (False, True)):
        """
        Constructor.

        Arguments:
            draw_size -- drawing area size,
            margins -- margins,
            pic_size -- picture size.
        """

        # Create image and drawer.
        self.Img = Image.new('RGB', pic_size, color = (230, 230, 230))
        self.Canvas = aggdraw.Draw(self.Img)
        self.Canvas.setantialias(True)

        # Data for transform.
        (dxl, dyl, dxh, dyh) = draw_area
        (px, py) = pic_size
        (mx, my) = margins
        self.FXI = (dxl, dxh)
        self.FYI = (dyl, dyh)
        (inv_x, inv_y) = invert
        if inv_x:
            self.TXI = (px - mx, mx)
        else:
            self.TXI = (mx, px - mx)
        if inv_y:
            self.TYI = (py - my, my)
        else:
            self.TYI = (my, py - my)

#---------------------------------------------------------------------------------------------------

    def FSS(self,
            filename = 'test.png'):
        """
        Flush, Save, Show.

        Arguments:
            filename -- name of file.
        """

        # Flush
        self.Canvas.flush()
        self.Img.save(filename)
        self.Img.show()

#---------------------------------------------------------------------------------------------------

    def TransformCoord(self, x, f, t):
        """
        Transform coordinate.

        Arguments:
            x -- coordinate,
            f -- from interval,
            t -- to interval.

        Result:
            New coordinate.
        """

        (fl, fh) = f
        (tl, th) = t
        a = (tl - th) / (fl - fh)
        b = tl - a * fl

        return a * x + b

#---------------------------------------------------------------------------------------------------

    def To(self, p):
        """
        Transform point TO.

        Arguments:
            p -- point from drawing area.

        Result:
            Point from picture area.
        """

        (px, py) = p

        return (self.TransformCoord(px, self.FXI, self.TXI),
                self.TransformCoord(py, self.FYI, self.TYI))

#---------------------------------------------------------------------------------------------------

    def From(self, p):
        """
        Transfrom point FROM.

        Arguments:
            p -- point from picture area.

        Result:
            Point from drawing area.
        """

        (px, py) = p

        return (self.TransformCoord(px, self.TXI, self.FXI),
                self.TransformCoord(py, self.TYI, self.FYI))

#---------------------------------------------------------------------------------------------------

    def Line(self, p1, p2,
             pen = aggdraw.Pen('black', 1.0)):
        """
        Line.

        Arguments:
            p1 -- from point,
            p2 -- to point,
            pen -- pen.
        """

        self.Canvas.line(self.To(p1) + self.To(p2), pen)

#---------------------------------------------------------------------------------------------------

    def FixLine(self, p, dp,
                pen = aggdraw.Pen('black', 1.0)):
        """
        Fix line.

        Arguments:
            p -- from point,
            dp -- fixed displacement.
        """

        (cx, cy) = self.To(p)
        (dpx, dpy) = dp
        c = (cx, cy, cx + dpx, cy + dpy)
        self.Canvas.line(c, pen)

#---------------------------------------------------------------------------------------------------

    def Ellipse(self, p1, p2,
                pen = aggdraw.Pen('black', 1.0),
                brush = None):
        """
        Ellipse.

        Arguments:
            p1 -- first point,
            p2 -- second point,
            pen -- pen,
            brush -- brush.
        """

        c = self.To(p1) + self.To(p2)
        if brush == None:
            self.Canvas.ellipse(c, pen)
        else:
            self.Canvas.ellipse(c, pen, brush)

#---------------------------------------------------------------------------------------------------

    def Point(self, p, r,
              pen = aggdraw.Pen('black', 1.0),
              brush = None):
        """
        Point.

        Arguments:
            p -- center of point,
            r -- radius,
            pen -- pen,
            brush -- brush.
        """

        (cx, cy) = self.To(p)
        c = (cx - r, cy - r, cx + r, cy + r)
        if brush == None:
            self.Canvas.ellipse(c, pen)
        else:
            self.Canvas.ellipse(c, pen, brush)

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

def test_drawer():
    """
    Drawer test.
    """

    D = Drawer()
    D.Line((0, 0), (100, 100))
    D.Ellipse((30, 30), (70, 70), pen = aggdraw.Pen('red', 1.5), brush = aggdraw.Brush('blue'))
    D.Point((10, 80), 5)
    D.Point((80, 10), 5)
    D.FSS()

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    #test_pil()
    #test_aggdraw()
    test_drawer()

#---------------------------------------------------------------------------------------------------
