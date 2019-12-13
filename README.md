# Raytrace3

To incorporate movement, we’ll need to a new file; a script file. This is independent of the scene file.
This file will only contain two types of entries: translation and rotation.
translation entries look like this:
0.000000x 0.000000y 5.972222z
This means: Move 0 in X, 0 in Y, 5.972222 in Z. Note that the letters x, y, and z in the output file are lowercase.
Or this:
0.000000rx 0.750000ry 0.000000rz 0XO 0YO 800ZO
This is a rotation. rx is a rotation about x, ry is a rotation about y, and rz is a rotation about z. Note that they are lower case. Also on this line is the point you want to shift the origin to. This uses XO (letters x and o), YO, and ZO to indicate the “new” origin. This particular line means treat 0,0,800 as the origin.
