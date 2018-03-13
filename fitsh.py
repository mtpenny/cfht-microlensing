import numpy as np

def eval_2d_poly(x,y,order,coeff,ox,oy,scale):

    order=int(order)
    x=(x-ox)/scale
    y=(y-oy)/scale
    if order==0:
        return coeff[0]
    elif order==1:
        return coeff[0]+x*coeff[1]+y*coeff[2]
    else:
        ret=0.0
        m=order*(order+1)/2
        n=order
        ym=1.0
        yf=1
        while n>=0:
            i=n
            p=m
            w=coeff[m]
            k=order
            while i>=1:
                w=w*x/i
                p=p-k
                w=w+coeff[p]
                i=i-1
                k=k-1
            ret=ret+w*ym
            ym=ym*y/yf
            yf=yf+1
            n=n-1
            m=m+1
        return ret

def get_jacobi(dxfit,dyfit,params):

    order = params[0]
    offsetx = params[1]
    offsety = params[2]
    scale = params[3]

    jnvar = order*(order+1)/2
    jxx = np.zeros(jnvar)
    jxy = np.zeros(jnvar)
    jyx = np.zeros(jnvar)
    jyy = np.zeros(jnvar)

    k=0
    for i in range(order):
        for j in range(i+1):
            jxx[k] = jxx[k] + dxfit[k+i+1]/scale
            jxy[k] = jxy[k] + dxfit[k+i+2]/scale
            jyx[k] = jyx[k] + dyfit[k+i+1]/scale
            jyy[k] = jyy[k] + dyfit[k+i+1]/scale
            k=k+1

    return [jxx,jxy,jyx,jyy]
        
def invertpoly2d(x,y,dxfit,dyfit,params):
    #Apply fitsh's inverse polynomial transform
    #fitsh function transformation_eval_invert_2d
    order = params[0]
    offsetx = params[1]
    offsety = params[2]
    scale = params[3]

    wx = x-dxfit[0]
    wy = y-dyfit[0]
    mxx = dxfit[1]
    mxy = dxfit[2]
    myx = dyfit[1]
    myy = dyfit[2]
    det = 1.0/(mxx*myy-mxy*myx)
    imxx = +myy*det
    imxy = -mxy*det
    imyx = -myx*det
    imyy = +mxx*det
    x0 = imxx*wx + imxy*wy
    y0 = imyx*wx + imyy*wy

    jxx,jxy,jyx,jyy = get_jacobi(dxfit,dyfit,params)

    n = 100
    px0 = x0
    py0 = y0
    if order>=2:
        for i in range(n):
            mxx = eval_2d_poly(x0,y0,order-1,jxx,offsetx,offsety,scale)
            mxy = eval_2d_poly(x0,y0,order-1,jxy,offsetx,offsety,scale)
            myx = eval_2d_poly(x0,y0,order-1,jyx,offsetx,offsety,scale)
            myy = eval_2d_poly(x0,y0,order-1,jyy,offsetx,offsety,scale)
            det=1.0/(mxx*myy-mxy*myx)
	    imxx=+myy*det
	    imxy=-mxy*det
	    imyx=-myx*det
	    imyy=+mxx*det
            wx=eval_2d_poly(x0,y0,order,dxfit,offsetx,offsety,scale)-x
            wy=eval_2d_poly(x0,y0,order,dyfit,offsetx,offsety,scale)-y
            dx=imxx*wx+imxy*wy
	    dy=imyx*wx+imyy*wy
	    x0=x0-dx
	    y0=y0-dy

    return x0,y0
    
def projection_get_matrix(ra0,dec0):
    mproj = np.zeros(shape=[3,3])
    ra=ra0*np.pi/180.0
    dec=dec0*np.pi/180.0

    sr0=np.sin(ra)
    cr0=np.cos(ra)
    sd0=np.sin(dec)
    cd0=np.cos(dec)

    mproj[0,0] = + sr0
    mproj[0,1]=-cr0
    mproj[0,2]=0.0
    mproj[1,0]=-sd0*cr0
    mproj[1,1]=-sd0*sr0
    mproj[1,2]=+cd0
    mproj[2,0]=-cd0*cr0
    mproj[2,1]=-cd0*sr0
    mproj[2,2]=-sd0

    return mproj

def projection_do_matrix_coord(mproj,rain,decin):
    ra = rain*np.pi/180.0
    dec = decin*np.pi/180.0

    sr=np.sin(ra)
    cr=np.cos(ra)
    sd=np.sin(dec)
    cd=np.cos(dec)

    x=cd*cr
    y=cd*sr
    z=sd

    rx=mproj[0,0]*x+mproj[0,1]*y+mproj[0,2]*z
    ry=mproj[1,0]*x+mproj[1,1]*y+mproj[1,2]*z
    rz=mproj[2,0]*x+mproj[2,1]*y+mproj[2,2]*z

    m=1.0/np.sqrt(1-rx**2-ry**2) #gnomic projection
    rx = rx*m
    ry = ry*m

    return [-rx*180.0/np.pi,ry*180.0/np.pi]

def projection_do_inverse_matrix_coord(mproj,x,y):

    z=1.0-x**2-y**2
    z=-np.sqrt(z);

    px=mproj[0,0]*x+mproj[1,0]*y+mproj[2,0]*z
    py=mproj[0,1]*x+mproj[1,1]*y+mproj[2,1]*z
    pz=mproj[0,2]*x+mproj[1,2]*y+mproj[2,2]*z

    rde=np.asin(pz)*180.0/np.pi
    rra=np.atan2(py,px)*180.0/np.pi
    if rra<0.0:
        rra=rra+360.0

    return [rra,rdec]

def ad2pix(ra,dec,astrom_data,dxfit,dyfit):

    params = astrom_data
    
    #Compute the projection matrix
    mproj = projection_get_matrix(params[-2],params[-1])

    #Project the sky coordinates onto the tangent plane
    rx,ry =  projection_do_matrix_coord(mproj,ra,dec)
    
    #Transform the projected coordinates to the pixel coordinates
    px = eval_2d_poly(rx,ry,params[0],dxfit,params[1],params[2],params[3])
    py = eval_2d_poly(rx,ry,params[0],dyfit,params[1],params[2],params[3])

    return np.array([px,py]).T
