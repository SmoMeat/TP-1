from turtle import *

def draw_square(x, y, side_length=15):
    """For visual studio"""
    penup(), goto(x, y), pendown()
    write('  '+'M', move=False, align='left', font=('Arial', 8, 'normal'))
    for side in range(4):
        fd(side_length), lt(90)
        
draw_square(0,0)

mainloop()