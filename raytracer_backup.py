"""
RayTracer Project : 
This project is to create rays that intersects with a sphere/lens and then 
refracts and intersects at an output plane. 
p is the starting point of the ray
k is direction of the ray
1 ray is created and then for loop is used to create a bundle of rays.
A plano convex lens is created and a bundle of rays is propagated through it.
A spherical reflecting surface was created.
Rays are propagated to allow reflection on its surface.

by Eunice Chen 21/3/2014
"""

import numpy as np
import matplotlib.pyplot as plt

class Ray:
    def __init__(self,p=[0.,0.,0.],k=[0.,0.,0.]):
        ''' This generates a point and its direction.
            p = point
            k = direction
        '''
        self._p=np.array(p)
        self._k=np.array(k)
        self._listofpointsp = [p]
        
    def p(self):
        return self._p
        
    def k(self):
        return self._k
        
    def append(self,new_p,new_k):
        '''
        creating a new point and a new direction for the ray
        '''
        self._p = new_p
        self._k = new_k
        return self._listofpointsp.append(new_p)

    def vertices(self):
        '''
        returning the old p and new p together
        '''
        return self._listofpointsp

''' This is a base class
'''
class OpticalElement:

    def propagate_ray(self, ray):
        '''
        It generates an optical ray at initial point, p and direction k
        that propagates through the optical element
        '''
        raise NotImplementedError()
    
''' This is the derived class of Optical Element. 
    This class represents a spherical surface.
    A ray is propagated through this class. Refraction occurs when the ray intersects the spherical surface.
'''
class SphericalRefraction(OpticalElement):
    
    def __init__(self,z0,curv,n1,n2,apradius):
        '''
        It generates an optical ray that refracts.
        z0 = intercept of the surface with the z-axis
        curv = curvature of the surface, 
        n1, n2 = refractive indices either side of the surface 
        apradius = maximum extent of the surface from the optical axis.
        '''
        self.z0=z0
        self.curv=curv
        self.n1=n1
        self.n2=n2
        self.apradius=apradius #This is aperture radius, it extends from the optical axis, in the direction perpendicular to the optical axis.
        self.R = 0 # making 0 a default value for self.R. When curv is 0, the R also becomes 0 to avoid an error

        try:
            self.R = 1/curv
        except ZeroDivisionError:
            print "R becomes 0"
        
        self.origin=np.array([0.0,0.0,self.z0 + self.R])
        
    def intercept(self,ray):
        '''
        This generates a ray that intersects with a spherical surface.
        It creates conditions of where the ray intersects twice, once or missed the sphere.
        The conditions for aperture radius is then created, where if x**2+y**2 >= aperture radius, the ray will miss the lens.
        '''
        r= ray.p() - self.origin
        magnitude_r=np.linalg.norm(r)
        k_hat = ray.k()/np.linalg.norm(ray.k())
        a=((np.dot(r,k_hat))**2)-((magnitude_r)**2-(self.R)**2)
        
        if a < 0:
            print "The ray missed the sphere"
            return None
            
        elif a == 0:
            l = -np.dot(r,k_hat)
            Q = ray.p()+(l*k_hat)
            return Q
            
        elif a > 0:
            l1=-np.dot(r,k_hat)+np.sqrt(a)
            l2=-np.dot(r,k_hat)-np.sqrt(a)
            Q1 = ray.p()+(l1*k_hat)
            Q2 = ray.p()+(l2*k_hat)
            
            if l1 < l2:
                
                if (Q1[0])**2+(Q1[1])**2 >= (self.apradius)**2:
                    print "The ray missed the lens - case 1"
                    return None
                    
                elif (Q1[0])**2+(Q1[1])**2 < (self.apradius)**2:
                    return Q1
            else:
                    
                if (Q2[0])**2+(Q2[1])**2 >= (self.apradius)**2:
                    print "The ray missed the lens - case 2"
                    return None
                    
                elif (Q2[0])**2+(Q2[1])**2 < (self.apradius)**2:
                    return Q2
                
        #The code should never reach this because all possibilities have been considered    
        else:
            print "Something unexpected happened"
            return None
        
    def refraction(self,ray):
        '''
        This generates a ray that refractes when it intersects the sphere.
        Total internal reflection is taken into account, so when that happens, it returns None
        The formula for snell's law refraction in vector form is used.
        '''
        ratio_n=self.n2/self.n1
        khati=ray.k()/np.linalg.norm(ray.k())
        surfnorm_n = self.intercept(ray) - self.origin
        surfnorm_normalised = surfnorm_n/np.linalg.norm(surfnorm_n)
        b=np.dot(surfnorm_normalised,khati)
        sinetheta1=np.sqrt(1-b**2)
        
        if sinetheta1 > ratio_n:
            print "Total internal reflection occured"
            return None
        
        c=np.dot(-surfnorm_normalised,khati)
        d=1/ratio_n
        e=1-((d**2)*(1-c**2))
        
        refracted_ray = (d*khati)+(((d*c)-np.sqrt(e))*surfnorm_normalised)
        return refracted_ray
        
    def propagate_ray(self,ray):
        '''
        This method propagates the ray so that it refracts when it intersects the sphere.
        When the intercept or refraction method returns None,
        the code is Terminated.
        Otherwise, a new point and direction, p1,k1 is created.
        '''
        if self.intercept(ray) is None or self.refraction(ray) is None:
            print "intercept or refraction is None"
            return "Terminated"
            
        else:
            p1=self.intercept(ray)
            k1=self.refraction(ray)
            ray.append(p1,k1)
            return ray.vertices()
'''
This class creates an infinite plane for the points to intersect.
'''
class OutputPlane(OpticalElement):
    
    def __init__(self,z):
        self.z=z
        
    def intercept(self,ray):
        '''
        This creates an intersection point on the output plane, 
        from the intersection point on the sphere.
        '''
        lamda=(self.z-ray.p()[2])/(ray.k()[2])
        y=ray.p()[1]+(lamda*ray.k()[1])
        x=ray.p()[0]+(lamda*ray.k()[0])
        OPintersection = ([x,y,self.z])
        return OPintersection
        
    def propagate_ray(self,ray):
        '''
        This propagates the ray from the sphere to the output plane.
        '''
        p2=self.intercept(ray)
        k2=ray.k()
        ray.append(p2, k2)
        return ray.vertices()
        
def trialplot(ray):
    '''
    This is to test if the code works.
    It plots one ray at a time.
    Plot two different rays that start at different points, if they converge after refraction,
    it shows that the code works.
    '''
    x=tuple(j[0] for j in ray.vertices()) #for each ray in the vectices take the x coordinate
    z=tuple(j[2] for j in ray.vertices()) #this takes the z coordinate
    plt.plot(z,x)

def rtpairs(R,N):
    '''generates r,n pairs
       r and n are lists of the same amount of numbers 
       r is the radii of the rings.
       R is a list of radii, r.
       N is a list containing the number of angles for each value or r.
    '''
    for i in range(len(R)):
	r=R[i]
	n=N[i]
        theta = 0.0 #starting with theta = 0
        for j in range(N[i]): 
            theta +=2*(np.pi)/N[i]
            yield R[i],theta
            
def rtuniform(n, rmax, m):
    '''generate r,n with uniform spacing
       n is the number of rings plotted
       rmax is the maximum radius of the plot
       m is the number of points each ring is increased by.
    '''
    R=np.arange(0,rmax,rmax/n)
    N=np.arange(1, n*m, m)
    return rtpairs(R, N)
    
def startingpoint(n,rmax,m):
    for r,t in rtuniform(n, rmax, m):
        yield np.array([r*np.cos(t), r*np.sin(t),0.])

'''
This class is to allow a refraction of a ray through a plane instead of a sphere
Most of the code is similar to the one in SphericalRefraction class,
but in the refraction method, the surface normal is different.
'''

class Planerefraction(OpticalElement):
    def __init__(self, z0, n2, n1, apradius):
        self.z0=z0
        self.n1=n1
        self.n2=n2
        self.apradius=apradius
        
    def intercept(self,ray):
        
        lamda=(self.z0-ray.p()[2])/(ray.k()[2])
        y=ray.p()[1]+(lamda*ray.k()[1])
        x=ray.p()[0]+(lamda*ray.k()[0])
        OPintersection = ([x,y,self.z0])
        return OPintersection
        
    def refraction(self,ray):
        '''
        This generates a ray that refractes when it intersects the plane.
        Total internal reflection is taken into account, so when that happens, it returns None
        The formula for snell's law refraction in vector form is used.
        '''
        vector1=np.array([1.,1.,self.z0])
        vector2=np.array([-1.,1.,self.z0])
        r=np.array([0.,0.,self.z0])
        dir_vec1=np.array([vector1[0]-r[0],vector1[1]-r[1],vector1[2]-r[2]])
        dir_vec2=np.array([vector2[0]-r[0],vector2[1]-r[1],vector2[2]-r[2]])
        ratio_n=self.n1/self.n2
        khati=ray.k()/np.linalg.norm(ray.k())
        surfnorm_n = -np.cross(dir_vec1,dir_vec2)
        surfnorm_normalised = surfnorm_n/np.linalg.norm(surfnorm_n)
        b=np.dot(surfnorm_normalised,khati)
        sinetheta1=np.sqrt(1-b**2)
        
        if sinetheta1 > ratio_n:
            print "Total internal reflection occured"
            return None
        
        c=np.dot(-surfnorm_normalised,khati)
        d=1/ratio_n
        e=1-((d**2)*(1-c**2))
        
        refracted_ray = (d*khati)+(((d*c)-np.sqrt(e))*surfnorm_normalised)
        return refracted_ray
        
    
    def propagate_ray(self,ray):
        '''
        This method propagates the ray so that it refracts when it intersects the plane.
        When the intercept or refraction method returns None,
        the code is Terminated.
        Otherwise, a new point and direction, p1,k1 is created.
        '''
        if self.refraction(ray) is None:
            print "refraction is None"
            return "Terminated"
            
        else:
            p2=self.intercept(ray)
            k2=self.refraction(ray)
            ray.append(p2,k2)
            return ray.vertices()

'''
This class creates a Plano Convex Lens. 
It takes the ray through the methods in the Spherical Refraction class
and then to the Planerefraction class. This makes the ray refract twice.
It refracts with a spherical surface first, then refracts with a plane surface and vice versa.
'''
        
class PlanoConvexLens:
    def __init__(self, z0, curv, n1, n2, apradius):

        self.z0=z0
        self.curv=curv
        self.n1=n1
        self.n2=n2
        self.apradius=apradius
        self.curvedsurface = SphericalRefraction(self.z0, self.curv, self.n1, self.n2, self.apradius)
        self.plane = Planerefraction(self.z0, self.n2, self.n1, self.apradius)
        
    def propagate_ray(self,ray):
        s1=self.curvedsurface.propagate_ray(ray) 
        return ray.vertices()
        
'''
This class outputs the refracted ray through the plano convex lens onto an output plane.
'''
class Output:
    def __init__(self,z):
        self.z=z
        self.output=OutputPlane(self.z)
        
    def intercept(self,ray):
        '''
        This creates an intersection point on the output plane, 
        from the intersection point on the sphere.
        '''
        lamda=(self.z-ray.p()[2])/(ray.k()[2])
        y=ray.p()[1]+(lamda*ray.k()[1])
        x=ray.p()[0]+(lamda*ray.k()[0])
        OPintersection = ([x,y,self.z])
        return OPintersection
    
    def propagate_ray(self,ray):
        s3=self.output.propagate_ray(ray)
        return ray.vertices()
    
'''
This class represents a spherical surface.
Reflection occurs when a ray is propagated through this class
'''
class SphericalReflection(OpticalElement):
    def __init__(self,z0,curv,apradius):
        self.z0=z0
        self.curv=curv
        self.apradius=apradius #This is aperture radius, it extends from the optical axis, in the direction perpendicular to the optical axis.
        self.R = 0 # making 0 a default value for self.R. When curv is 0, the R also becomes 0 to avoid an error

        try:
            self.R = 1/curv
        except ZeroDivisionError:
            print "R becomes 0"
        
        self.origin=np.array([0.0,0.0,self.z0 + self.R])
    
    def intercept(self,ray):
        '''
        This generates a ray that intersects with a spherical surface.
        It creates conditions of where the ray intersects twice, once or missed the sphere.
        The conditions for aperture radius is then created, where if x**2+y**2 > (aperture radius)**2, the ray will miss the lens.
        '''
        r= ray.p() - self.origin
        magnitude_r=np.linalg.norm(r)
        k_hat = ray.k()/np.linalg.norm(ray.k())
        a=((np.dot(r,k_hat))**2)-((magnitude_r)**2-(self.R)**2)
        
        if a < 0:
            print "The ray missed the sphere"
            return None
            
        elif a == 0:
            l = -np.dot(r,k_hat)
            Q = ray.p()+(l*k_hat)
            print 'There is 1 intersection'
            return Q
            
        elif a > 0:
            l1=-np.dot(r,k_hat)+np.sqrt(a)
            l2=-np.dot(r,k_hat)-np.sqrt(a)
            Q1 = ray.p()+(l1*k_hat)
            Q2 = ray.p()+(l2*k_hat)
            
            if l1 < l2:
                
                if (Q1[0])**2+(Q1[1])**2 > (self.apradius)**2:
                    print "The ray missed the lens - case 1"
                    return None
                    
                elif (Q1[0])**2+(Q1[1])**2 <= (self.apradius)**2:
                    return Q1
            else:
                    
                if (Q2[0])**2+(Q2[1])**2 > (self.apradius)**2:
                    print "The ray missed the lens - case 2"
                    return None
                    
                elif (Q2[0])**2+(Q2[1])**2 <= (self.apradius)**2:
                    return Q2
                
        #The code should never reach this because all possibilities have been considered    
        else:
            print "Something unexpected happened"
            return None
            
    def reflection(self,ray):
        
        khati=ray.k()/np.linalg.norm(ray.k())
        surfnorm_n = self.intercept(ray) - self.origin
        surfnorm_normalised = surfnorm_n/np.linalg.norm(surfnorm_n)
        reflected_ray = khati + (2*(np.dot(-surfnorm_normalised,khati))*surfnorm_normalised)
        return reflected_ray
        
    def propagate_ray(self,ray):
        '''
        This method propagates the ray so that it reflects when it intersects the sphere.
        When the intercept or reflection method returns None,
        the code is Terminated.
        Otherwise, a new point and direction, p1,k1 is created.
        '''
        if self.intercept(ray) is None or self.reflection(ray) is None:
            print "intercept or reflection is None"
            return "Terminated"
            
        else:
            p1=self.intercept(ray)
            k1=self.reflection(ray)
            ray.append(p1,k1)
            return ray.vertices()
