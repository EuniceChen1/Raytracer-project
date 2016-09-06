For ray bundle with normal propagation : #do everything in ipython
for i in arange(-1.,1.,0.05):
    a=raytracer.OutputPlane(250)
    b=raytracer.Ray(p=[i,0.,0.],k=[0.,0.,1.])
    c=raytracer.SphericalRefraction(100.,0.03,1.,1.5,10.)
    return_value = c.propagate_ray(b)
    if return_value == "Terminated":
        print "One ray was terminated"
    else:
        a.propagate_ray(b)
        x=(b.vertices()[0][0],b.vertices()[1][0],b.vertices()[2][0])
        z=(b.vertices()[0][2],b.vertices()[1][2],b.vertices()[2][2])
        plt.plot(z,x)
        plt.title('Z-X plane')
        plt.xlabel('z')
        plt.ylabel('x')
        plt.show()
		
Another example for normal propagation ray bundle:
for i in arange(-3.,3.,0.3):
    a=raytracer.OutputPlane(200)
    b=raytracer.Ray(p=[i,0.,0.],k=[0.,0.,1.])
    c=raytracer.SphericalRefraction(100.,0.03,1.,1.5,10.)
    return_value = c.propagate_ray(b)
    if return_value == "Terminated":
        print "One ray was terminated"
    else:
        a.propagate_ray(b)
        x=(b.vertices()[0][0],b.vertices()[1][0],b.vertices()[2][0])
        z=(b.vertices()[0][2],b.vertices()[1][2],b.vertices()[2][2])
        plt.plot(z,x,'b-')
        plt.title('Z-X plane')
        plt.xlabel('z')
        plt.ylabel('x')
        plt.show()
		
For ray bundle propagated in an angle:
for i in arange(-3.,3.,0.3):
    a=raytracer.OutputPlane(200)
    b=raytracer.Ray(p=[i+10.,0.,0.],k=[-13.,-5.,100.])
    c=raytracer.SphericalRefraction(100.,0.03,1.,1.5,10.)
    return_value = c.propagate_ray(b)
    if return_value == "Terminated":
        print "One ray was terminated"
    else:
        a.propagate_ray(b)
        x=(b.vertices()[0][0],b.vertices()[1][0],b.vertices()[2][0])
        z=(b.vertices()[0][2],b.vertices()[1][2],b.vertices()[2][2])
        plt.plot(z,x,'r-')
        plt.title('Z-X plane')
        plt.xlabel('z')
        plt.ylabel('x')
        plt.show()


		
For initial ray with normal propagation:
for start_point in raytracer.startingpoint(10,10.,6):
    x0 = start_point[0]
    y0 = start_point[1]
    a=raytracer.OutputPlane(250)
    b=raytracer.Ray(p=[x0, y0, 0.],k=[0.,0.,1.])
    c=raytracer.SphericalRefraction(100.,0.03,1.,1.5,10.)
    return_value = c.propagate_ray(b)
    if return_value == "Terminated":
        print "One ray was terminated"
    else:
        intersectpoint=a.intercept(b)
        a.propagate_ray(b)
        x=intersectpoint[0]
        y=intersectpoint[1]
        plot(x,y,'b.')
        plt.title('X-Y plane')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
		
For initial ray with angled propagation:
for start_point in raytracer.startingpoint(10,5.,4):
    x0 = start_point[0]+2.
    y0 = start_point[1]+2.
    a=raytracer.OutputPlane(250)
    b=raytracer.Ray(p=[x0, y0, 0.],k=[-8,0.,200.])
    c=raytracer.SphericalRefraction(100.,0.03,1.,1.5,10.)
    return_value = c.propagate_ray(b)
    if return_value == "Terminated":
        print "One ray was terminated"
    else:
        intersectpoint=a.intercept(b)
        a.propagate_ray(b)
        x=intersectpoint[0]
        y=intersectpoint[1]
        plot(x,y,'r.')
        plt.title('X-Y plane')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
		
For output ray with normal propagation:
for start_point in raytracer.startingpoint(10,10.,6):
    x0 = start_point[0]
    y0 = start_point[1]
    a=raytracer.OutputPlane(200)
    b=raytracer.Ray(p=[x0, y0, 0.],k=[0.,0.,1.])
    c=raytracer.SphericalRefraction(100.,0.03,1.,1.5,10.)
    return_value = c.propagate_ray(b)
    if return_value == "Terminated":
        print "One ray was terminated"
    else:
        intersectpoint=a.intercept(b)
        a.propagate_ray(b)
        x=intersectpoint[0]
        y=intersectpoint[1]
        plot(x,y,'b.')
        plt.title('X-Y plane')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
		
For output ray with angled propagation: (some rays missed the lens)
for start_point in raytracer.startingpoint(10,8.,12):
    x0 = start_point[0]+3.
    y0 = start_point[1]+0.
    a=raytracer.OutputPlane(197.5)
    b=raytracer.Ray(p=[x0, y0, 0.],k=[-8,0.,113.])
    c=raytracer.SphericalRefraction(100.,0.03,1.,1.5,10.)
    return_value = c.propagate_ray(b)
    if return_value == "Terminated":
        print "One ray was terminated"
    else:
        intersectpoint=a.intercept(b)
        a.propagate_ray(b)
        x=intersectpoint[0]
        y=intersectpoint[1]
        plot(x,y,'b.')
        plt.title('X-Y plane')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
		



