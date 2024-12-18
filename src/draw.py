"""
Drawing examples.

Created on Tue May  7 09:58:45 2019

@author: Rybakov
"""

import aggdraw
from PIL import Image, ImageDraw, ImageFont

#---------------------------------------------------------------------------------------------------
# Variables.
#---------------------------------------------------------------------------------------------------

BlackPen = aggdraw.Pen('black', 1.0)

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
                 invert = (False, True),
                 backcolor = (255, 255, 255)):
        """
        Constructor.

        Arguments:
            draw_size -- drawing area size,
            margins -- margins,
            pic_size -- picture size.
        """

        # Create image and drawer.
        self.Img = Image.new('RGB', pic_size, color = backcolor)
        self.Backcolor = backcolor
        self.Canvas = aggdraw.Draw(self.Img)
        self.Canvas.setantialias(True)

        # Data for transform.
        self.DrawArea = draw_area
        (dxl, dyl, dxh, dyh) = draw_area
        (px, py) = pic_size
        (mx, my) = margins
        self.FXI = (dxl, dxh)
        self.FYI = (dyl, dyh)
        self.TXI = (mx, px - mx)
        self.TYI = (my, py - my)
        (inv_x, inv_y) = invert
        self.InvKX = -1.0 if inv_x else 1.0
        self.InvKY = -1.0 if inv_y else 1.0

#---------------------------------------------------------------------------------------------------

    def fss(self,
            filename = 'test.png'):
        """
        Flush, Save, Show.

        Arguments:
            filename -- name of file.
        """

        # Flush
        self.Canvas.flush()
        if not (filename is None):
            self.Img.save(filename)
        self.Img.show()

#---------------------------------------------------------------------------------------------------

    @staticmethod
    def transform_coord(x, f, t, inv_k):
        """
        Transform coordinate.

        Arguments:
            x -- coordinate,
            f -- from interval,
            t -- to interval,
            inv_k -- invert coefficient.

        Result:
            New coordinate.
        """

        (fl, fh) = f
        (tl, th) = t

        if inv_k > 0.0:
            a = (tl - th) / (fl - fh)
            b = tl - a * fl
        else:
            a = (th - tl) / (fl - fh)
            b = tl - a * fh

        return a * x + b

#---------------------------------------------------------------------------------------------------

    def to(self, p):
        """
        Transform point TO.

        Arguments:
            p -- point from drawing area.

        Result:
            Point from picture area.
        """

        (px, py) = p

        return (self.transform_coord(px, self.FXI, self.TXI, self.InvKX),
                self.transform_coord(py, self.FYI, self.TYI, self.InvKY))

#---------------------------------------------------------------------------------------------------

    def fr(self, p):
        """
        Transform point FROM.

        Arguments:
            p -- point from picture area.

        Result:
            Point from drawing area.
        """

        (px, py) = p

        return (self.transform_coord(px, self.TXI, self.FXI, self.InvKX),
                self.transform_coord(py, self.TYI, self.FYI, self.InvKY))

#---------------------------------------------------------------------------------------------------

    def line(self, p1, p2, pen = BlackPen):
        """
        Line.

        Arguments:
            p1 -- from point,
            p2 -- to point,
            pen -- pen.
        """

        self.Canvas.line(self.to(p1) + self.to(p2), pen)

#---------------------------------------------------------------------------------------------------

    def fix_line(self, p, dp, pen=BlackPen):
        """
        Fix line.

        Arguments:
            p -- from point,
            dp -- fixed displacement.
        """

        (cx, cy) = self.to(p)
        (dpx, dpy) = dp
        c = (cx, cy, cx + dpx, cy + dpy)
        self.Canvas.line(c, pen)

#---------------------------------------------------------------------------------------------------

    def ellipse(self, p1, p2, pen = BlackPen, brush = None):
        """
        Ellipse.

        Arguments:
            p1 -- first point,
            p2 -- second point,
            pen -- pen,
            brush -- brush.
        """

        c = self.to(p1) + self.to(p2)
        if brush is None:
            self.Canvas.ellipse(c, pen)
        else:
            self.Canvas.ellipse(c, pen, brush)

#---------------------------------------------------------------------------------------------------

    def point(self, p, r, pen = BlackPen, brush = None):
        """
        Point.

        Arguments:
            p -- center of point,
            r -- radius,
            pen -- pen,
            brush -- brush.
        """

        (cx, cy) = self.to(p)
        c = (cx - r, cy - r, cx + r, cy + r)
        if brush is None:
            self.Canvas.ellipse(c, pen)
        else:
            self.Canvas.ellipse(c, pen, brush)

#---------------------------------------------------------------------------------------------------

    def axis(self,
             h = 12, w = 4,
             pen = aggdraw.Pen('silver', 2.0)):
        """
        Draw axis.

        Arguments:
            h -- arrow height,
            w -- arrow width,
            pen -- pen.
        """

        (min_x, min_y, max_x, max_y) = self.DrawArea

        # OX
        if (min_y <= 0) and (max_y >= 0):
            self.line((min_x, 0), (max_x, 0), pen = pen)
            self.fix_line((max_x, 0), (-h * self.InvKX, w), pen = pen)
            self.fix_line((max_x, 0), (-h * self.InvKX, -w), pen = pen)

        # OY
        if (min_x <= 0) and (max_x >= 0):
            self.line((0, min_y), (0, max_y), pen = pen)
            self.fix_line((0, max_y), (w, -h * self.InvKY), pen = pen)
            self.fix_line((0, max_y), (-w, -h * self.InvKY), pen = pen)

#---------------------------------------------------------------------------------------------------

    def grid(self,
             deltas,
             pen = aggdraw.Pen('silver', 1.0)):
        """
        Draw grid.

        Arguments:
            deltas -- distances between adjacent lines,
            pen -- pen.
        """

        (min_x, min_y, max_x, max_y) = self.DrawArea

        (gx, gy) = deltas
        gx_cur = gx
        while gx_cur <= max_x:
            if gx_cur >= min_x:
                self.line((gx_cur, min_y), (gx_cur, max_y), pen = pen)
            gx_cur = gx_cur + gx
        gx_cur = -gx
        while gx_cur >= min_x:
            if gx_cur <= max_x:
                self.line((gx_cur, min_y), (gx_cur, max_y), pen = pen)
            gx_cur = gx_cur - gx
        gy_cur = gy
        while gy_cur <= max_y:
            if gy_cur >= min_y:
                self.line((min_x, gy_cur), (max_x, gy_cur), pen = pen)
            gy_cur = gy_cur + gy
        gy_cur = -gy
        while gy_cur >= min_y:
            if gy_cur <= max_y:
                self.line((min_x, gy_cur), (max_x, gy_cur), pen = pen)
            gy_cur = gy_cur - gy

#---------------------------------------------------------------------------------------------------

    def rect(self, p1, p2, pen = BlackPen, brush = None):
        """
        Rectangle.

        Arguments:
            p1 -- first point,
            p2 -- second point,
            pen -- pen,
            brish -- brush.
        """

        if brush is None:
            self.Canvas.rectangle(self.to(p1) + self.to(p2), pen)
        else:
            self.Canvas.rectangle(self.to(p1) + self.to(p2), pen, brush)

#---------------------------------------------------------------------------------------------------

    def rect_with_center_point(self, p1, p2, r, pen = BlackPen, brush = None):
        """
        Rectangle with point in its center.

        Arguments:
            p1 -- first point,
            p2 -- second point,
            r -- point radius,
            pen -- pen,
            brush -- brush.
        """

        (p1x, p1y) = p1
        (p2x, p2y) = p2

        self.rect(p1, p2, pen, brush)
        self.point((0.5 * (p1x + p2x), 0.5 * (p1y + p2y)), r, pen, brush)

#---------------------------------------------------------------------------------------------------

    def full_graph(self, ps, pen = BlackPen):
        """
        Draw full graph.

        Arguments:
            ps -- point list,
            pen -- pen.
        """

        n = len(ps)
        for i in range(n):
            for j in range(i + 1, n):
                self.line(ps[i], ps[j], pen)

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

    d = Drawer()
    d.line((0, 0), (100, 100))
    d.ellipse((30, 30), (70, 70), pen = aggdraw.Pen('red', 1.5), brush = aggdraw.Brush('steelblue'))
    d.point((10, 80), 5)
    d.point((80, 10), 5)
    d.fss()

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test_pil()
    test_aggdraw()
    test_drawer()

#---------------------------------------------------------------------------------------------------
