Circlipse v.1.0, for Adobe Illustrator
======================================

**Circlipse** is an Illustrator script that enables quick and precise
drawing of ellipses and circles, even when a traced ellipse is partly
obscured. It is inspired by
[a cool Inkscape plugin](http://pernsteiner.org/inkscape/ellipse_5pts/)
and is intended for
[similar uses](http://www.youtube.com/watch?v=NAl3WJBT8Z8).
However, it is a purely custom implementation in Javascript
that emphasizes precision. It is heavily commented for anyone interested
in the source code.

How to Use
----------

1. Copy the JSX script file to the scripts folder. If you can find your
   Illustrator folder, go from there to `Presets` and then to `Scripts`.
   (The `Scripts` folder may itself be wrapped under a locale-named folder,
   like `en_US`.)
2. Start (or restart) Illustrator and open the target document. Save
   any work, just in case!
3. Set the stroke and fill properties you want.
4. With the Pen Tool, plot either three anchor points for a circle or five
   for an ellipse, with each point lying at the outline of the path desired.
5. Select File -> Scripts -> circlipse. If the script calculates the
   ellipse or circle correctly, you're done! Otherwise, pay close
   attention to any alert messages.

Further Development
-------------------

This script may be useful for other Adobe products, so if requested, I can
attempt making it run elsewhere.

Recent Updates
--------------
An older version of this script was prone to failing if the top layer were
hidden. I do not recall it being documented anywhere, but the active layer
has properties in the API similar to that of the document itself. Thus, the
current layer is now manipulated directly in lieu of adding an ellipse to
the top layer as was happening before, which should improve workflow. I have
also added the MIT license and confirmed on my copy of CS6 a while back.
