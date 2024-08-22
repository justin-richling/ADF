class MyClass:
    def __init__(self, my_var):
        self.my_var = my_var
    
    def outer_function(self):
        self.inner_function()
    
    def inner_function(self):
        # Access the variable via the object's state
        print(self.my_var)

obj = MyClass("some value")
obj.outer_function()
