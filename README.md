# A RayTracing render application
### Support box, cylinder, cone and sphere objects.
### Each object could have its own material parameter, reflection and transparency parameter.
### A scene could have multiple lights. Each light have its own position , direction and color. Natural light and spot light are supported.
### Scene is constructed by scene graph. Position, parameters of objects can be written in config file and animation code can be easily added to scene

## How to Run
You can run our program on different scenegraph models 
To change scenegraph models, you can change the input XML file inside JOGLFrame.java line 53. 
To capture ray tracing images, you can press "c" and the output image will be named raytracing.jpg and stored on D drive. When the program start running, it will also capture the first frame and save it. 
Each time you press "c", the program will capture a new image with the same name and stored on the same location. It means you press "c", raytracing.jpg is different. 
 
*based on utilites provided by prof. Amit Shesh@Northeastern University
