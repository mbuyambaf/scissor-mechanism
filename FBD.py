from logging import root
from tkinter import *
from PIL import ImageTk,Image

def slider_changed(event):
    #print(horizontal.get())
    return horizontal.get()

root = Tk()
canvas  = Canvas()
base  = Canvas()

root.title("Force diagram")
root.geometry("400x400")

#table
height = slider_changed
print(height)
table = canvas.create_line(50, height, 300, height, width=5)
canvas.pack()
canvas.move(table, 50, height)

horizontal = Scale(root, from_=30, to=200,command=slider_changed, orient=HORIZONTAL)
horizontal.pack()




#canvas.pack()

#base
#base.create_line(50, 200, 300, 200, width=5)
#base.pack()

root.mainloop()

