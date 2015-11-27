import tkinter as tk
import ipdb


class MainGui(tk.Frame):

    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.parent = parent
        self.pack()

        self.parent.geometry('450x450+500+500')
        self.parent.title('Welcome to repESP !')

        w = tk.Label(self, text='Hello, world!')
        w.pack()


if __name__ == '__main__':

    tk_root = tk.Tk()
    MainGui(tk_root)
    tk_root.mainloop()
