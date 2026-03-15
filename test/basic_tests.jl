A = Operator("X11")
B = Operator("Y11")
C = im * Operator("Z11")
@assert norm(A * B - C) == 0

A = Operator("XX1")
B = Operator("Y11")
C = im * Operator("ZX1")
@assert norm(A * B - C) == 0

A = Operator("XYZ") + Operator("111") + Operator("XXX")
B = Operator("YY1")
C = Operator("YY1") - Operator("ZZX") + im * Operator("Z1Z")
@assert norm(A * B - C) == 0

A = Operator("XYZ") + Operator("111") + Operator("XXX")
B = Operator("111")
C = Operator("XYZ") + Operator("111") + Operator("XXX")
@assert norm(A * B - C) == 0

A = Operator("XYZ") + Operator("YXY") + Operator("YYY")
B = Operator("Y11") + Operator("X1X")
C = Operator("1XY") - Operator("ZXZ") + (1 + im) * Operator("1YY") + (-1 + im) * Operator("ZYZ")
@assert norm(A * B - C) == 0
