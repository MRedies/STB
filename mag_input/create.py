import numpy as np 
from scipy.interpolate import griddata
import sys

def read_file(fname):
    with open(fname) as f:
        _, x,y,z = f.readline().split(" ")
        x = int(x)
        y = int(y)
        z = int(z)

        data = np.zeros((x*y*z, 6))
        for i in range(x*y*z):
            data[i,:] = [float(num) for num in f.readline().split()]

        _, num_vec = f.readline().split()
        num_vec = int(num_vec)

        vec = []
        for i in range(num_vec):
            vec.append(f.readline())

        return [x,y,z], data, num_vec, vec

def cut_data(dim, data, z_val):
    z   = data[:,2]
    sel = np.where(z>z_val)

    data   = data[sel,:][0,:,:]
    z      = data[:,2]
    dim[2] = np.unique(z).shape[0]
       
    return dim, data

def stretch(dim, data, factor):
    pos = data[:,:3]
    m_x = data[:,3]
    m_y = data[:,4]
    m_z = data[:,5]

    z_max = np.max(pos[:,2])

    new_space = np.zeros(pos.shape)
    print(new_space.shape)
    new_space[:,0] = pos[:,0] / factor
    new_space[:,1] = pos[:,1] / factor     
    new_space[:,2] = z_max - (z_max - pos[:,2] )/factor
    
    interp_m_x = griddata(pos, m_x, new_space, method="linear")
    interp_m_y = griddata(pos, m_y, new_space, method="linear")
    interp_m_z = griddata(pos, m_z, new_space, method="linear")

    new_dat = np.zeros(data.shape)

    new_dat[:,:3] = pos
    new_dat[:,3]  = interp_m_x
    new_dat[:,4]  = interp_m_y
    new_dat[:,5]  = interp_m_z

    for i in range(np.prod(dim)):
        new_dat[i,3:6] /= np.linalg.norm(new_dat[i,3:6])

    return new_dat

def shrink(dim, data, factor):
    print("shrink by {}".format(factor))
    pos = data[:,:3]
    m_x = data[:,3]
    m_y = data[:,4]
    m_z = data[:,5]

    z_max = np.max(pos[:,2])

    new_space = np.zeros(pos.shape)
    print(new_space.shape)
    new_space[:,0] = pos[:,0] * factor
    new_space[:,1] = pos[:,1] * factor     
    new_space[:,2] = z_max - (z_max - pos[:,2] ) * factor
    
    interp_m_x = griddata(pos, m_x, new_space, method="linear")
    interp_m_y = griddata(pos, m_y, new_space, method="linear")
    interp_m_z = griddata(pos, m_z, new_space, method="linear")

    oo_box = np.where(np.isnan(interp_m_x))
    interp_m_x[oo_box] = 0.0
    interp_m_y[oo_box] = 0.0
    interp_m_z[oo_box] = 1.0
    
    new_dat = np.zeros(data.shape)

    new_dat[:,:3] = pos
    new_dat[:,3]  = interp_m_x
    new_dat[:,4]  = interp_m_y
    new_dat[:,5]  = interp_m_z

    for i in range(np.prod(dim)):
        new_dat[i,3:6] /= np.linalg.norm(new_dat[i,3:6])

    return new_dat

def stretch_z(dim, data, factor):
    pos = data[:,:3]
    m_x = data[:,3]
    m_y = data[:,4]
    m_z = data[:,5]

    z_max = np.max(pos[:,2])

    new_space = np.zeros(pos.shape)
    print(new_space.shape)
    new_space[:,0] = pos[:,0] 
    new_space[:,1] = pos[:,1]     
    new_space[:,2] = z_max - (z_max - pos[:,2] )/factor
    
    interp_m_x = griddata(pos, m_x, new_space, method="linear")
    interp_m_y = griddata(pos, m_y, new_space, method="linear")
    interp_m_z = griddata(pos, m_z, new_space, method="linear")

    new_dat = np.zeros(data.shape)

    new_dat[:,:3]  = pos
    new_dat[:,3]   = interp_m_x
    sel            = np.where(np.isnan(new_dat[:,3]))
    new_dat[sel,3] = 0.0

    new_dat[:,4]   = interp_m_y
    sel            = np.where(np.isnan(new_dat[:,4]))
    new_dat[sel,4] = 0.0

    new_dat[:,5]   = interp_m_z
    sel            = np.where(np.isnan(new_dat[:,5]))
    new_dat[sel,5] = 1.0

    for i in range(np.prod(dim)):
        new_dat[i,3:6] /= np.linalg.norm(new_dat[i,3:6])




    return new_dat


def ferro_below(data, crit_z):
    z   = data[:,2]
    sel = np.where(z <= crit_z)

    data[sel,3:5] = 0.0
    data[sel,5]   = 1.0 

    return data

def cart_to_sphere(x,y,z):
    theta = np.arccos(z / np.sqrt(x**2 + y**2 + z**2))
    phi   = np.arctan2(y,x)
    
    return theta, phi

def sphere_to_cart(theta, phi):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return x,y,z
    
def mixed_state(alpha):
    dim, bobb_data, num_vec, vec  = read_file("bobber_large.txt")
    _,   skym_data, _,       _    = read_file("skyrm_tube.txt")
    skyrm_data = ferro_below(skym_data, 24)

    bobb_theta, bob_phi = cart_to_sphere(bobb_data[:,3],
                                         bobb_data[:,4],
                                         bobb_data[:,5],)
    skyrm_theta, skyrm_phi = cart_to_sphere(skyrm_data[:,3],
                                            skyrm_data[:,4],
                                            skyrm_data[:,5],)

    theta = alpha * bobb_theta + (1-alpha) * skyrm_theta
    phi   = alpha * bob_phi    + (1-alpha) * skyrm_phi

    x,y,z = sphere_to_cart(theta, phi)

    new_data       = np.zeros(bobb_data.shape)
    new_data[:,:3] = bobb_data[:,:3]
    new_data[:,3]  = x
    new_data[:,4]  = y
    new_data[:,5]  = z

    return dim, new_data, num_vec, vec

def make_globuli(center_z, bobb_file="bobber_large.txt"):
    dim, data, num_vec, vec  = read_file(bobb_file)

    x   = data[:,0].reshape(dim, order="F")
    y   = data[:,1].reshape(dim, order="F")
    z   = data[:,2].reshape(dim, order="F")
    m_x = data[:,3].reshape(dim, order="F")
    m_y = data[:,4].reshape(dim, order="F")
    m_z = data[:,5].reshape(dim, order="F")

    new_m_x = np.zeros(m_x.shape)
    new_m_y = np.zeros(m_y.shape)
    new_m_z = np.ones(m_z.shape)

    if(center_z is None):
        center_z = np.mean(z)
    z_idx = np.argmin(np.abs(center_z-z[0,0,:]))
    z_sz  = np.unique(z).shape[0]

    
    new_m_x[:,:,:z_idx] = m_x[:,:,-z_idx:]
    new_m_y[:,:,:z_idx] = m_y[:,:,-z_idx:]
    new_m_z[:,:,:z_idx] = m_z[:,:,-z_idx:]
    new_m_x[:,:,z_idx:] = np.flip(m_x[:,:,-(z_sz-z_idx):], axis=2)
    new_m_y[:,:,z_idx:] = np.flip(m_y[:,:,-(z_sz-z_idx):], axis=2)
    new_m_z[:,:,z_idx:] = np.flip(m_z[:,:,-(z_sz-z_idx):], axis=2)

    new_data = np.zeros(data.shape)
    new_data[:,0] = x.flatten(order="F")
    new_data[:,1] = y.flatten(order="F")
    new_data[:,2] = z.flatten(order="F")
    new_data[:,3] = new_m_x.flatten(order="F")
    new_data[:,4] = new_m_y.flatten(order="F")
    new_data[:,5] = new_m_z.flatten(order="F")

    return dim, new_data, num_vec, vec 


def dbl_bobber(cut_h=8):
    dim, data, num_vec, vec  = read_file("bobber_large.txt")
    dim, data = cut_data(dim, data, cut_h)

    x   = data[:,0].reshape(dim, order="F")
    y   = data[:,1].reshape(dim, order="F")
    z   = data[:,2].reshape(dim, order="F")
    m_x = data[:,3].reshape(dim, order="F")
    m_y = data[:,4].reshape(dim, order="F")
    m_z = data[:,5].reshape(dim, order="F")

    new_m_x = np.zeros(m_x.shape)
    new_m_y = np.zeros(m_y.shape)
    new_m_z = np.zeros(m_z.shape)

    center_z = np.mean(np.unique(z))
    print("center = {}".format(center_z))
    z_idx = np.argmin(np.abs(center_z-z[0,0,:]))
    z_sz  = np.unique(z).shape[0]
    
    new_m_x[:,:,:z_idx] = np.flip(m_x[:,:,-z_idx:], axis=2)
    new_m_y[:,:,:z_idx] = np.flip(m_y[:,:,-z_idx:], axis=2)
    new_m_z[:,:,:z_idx] = np.flip(m_z[:,:,-z_idx:], axis=2)
    new_m_x[:,:,z_idx:] = m_x[:,:,-(z_sz-z_idx):]
    new_m_y[:,:,z_idx:] = m_y[:,:,-(z_sz-z_idx):]
    new_m_z[:,:,z_idx:] = m_z[:,:,-(z_sz-z_idx):]

    new_data = np.zeros(data.shape)
    new_data[:,0] = x.flatten(order="F")
    new_data[:,1] = y.flatten(order="F")
    new_data[:,2] = z.flatten(order="F")
    new_data[:,3] = new_m_x.flatten(order="F")
    new_data[:,4] = new_m_y.flatten(order="F")
    new_data[:,5] = new_m_z.flatten(order="F")

    return dim, new_data, num_vec, vec 

def skyrm_layer(from_top, cut_h=17):
    dim, data, num_vec, vec  = read_file("skyrm_tube.txt")
    dim, data = cut_data(dim, data, cut_h)

    x   = data[:,0].reshape(dim, order="F")
    y   = data[:,1].reshape(dim, order="F")
    z   = data[:,2].reshape(dim, order="F")
    m_x = data[:,3].reshape(dim, order="F")
    m_y = data[:,4].reshape(dim, order="F")
    m_z = data[:,5].reshape(dim, order="F")

    new_m_x = np.zeros(m_x.shape)
    new_m_y = np.zeros(m_y.shape)
    new_m_z = np.zeros(m_z.shape)
    
    from_top += 1
    
    new_m_x[:,:,:] = 0.0
    new_m_y[:,:,:] = 0.0
    new_m_z[:,:,:] = 1.0
    new_m_x[:,:,-from_top] = m_x[:,:,-1]
    new_m_y[:,:,-from_top] = m_y[:,:,-1]
    new_m_z[:,:,-from_top] = m_z[:,:,-1]

    new_data = np.zeros(data.shape)
    new_data[:,0] = x.flatten(order="F")
    new_data[:,1] = y.flatten(order="F")
    new_data[:,2] = z.flatten(order="F")
    new_data[:,3] = new_m_x.flatten(order="F")
    new_data[:,4] = new_m_y.flatten(order="F")
    new_data[:,5] = new_m_z.flatten(order="F")

    return dim, new_data, num_vec, vec 

def mix_glob_tube(a):
    dim, glob_data, num_vec, vec = make_globuli(None, bobb_file="bobber_gt8.txt")

    tube_dim, tube_data, _,_  = read_file("skyrm_tube.txt")
    _, tube_data = cut_data(tube_dim, tube_data, 8)

   
    glob_theta, glob_phi = cart_to_sphere(glob_data[:,3],
                                         glob_data[:,4],
                                         glob_data[:,5],)
    tube_theta, tube_phi = cart_to_sphere(tube_data[:,3],
                                            tube_data[:,4],
                                            tube_data[:,5],)

    theta = a * glob_theta + (1-a) * tube_theta
    phi   = a * glob_phi    + (1-a) * tube_phi

    x,y,z = sphere_to_cart(theta, phi)

    new_data       = np.zeros(glob_data.shape)
    new_data[:,:3] = glob_data[:,:3]
    new_data[:,3]  = x
    new_data[:,4]  = y
    new_data[:,5]  = z

    return dim, new_data, num_vec, vec


def mix_dbl_tube(a):
    dim, dbl_data, num_vec, vec = read_file("dbl_bobber.txt")

    tube_dim, tube_data, _,_  = read_file("skyrm_tube.txt")
    _, tube_data = cut_data(tube_dim, tube_data, 8)

   
    dbl_theta, dbl_phi = cart_to_sphere(dbl_data[:,3],
                                         dbl_data[:,4],
                                         dbl_data[:,5],)
    tube_theta, tube_phi = cart_to_sphere(tube_data[:,3],
                                            tube_data[:,4],
                                            tube_data[:,5],)

    theta = a * dbl_theta + (1-a) * tube_theta
    phi   = a * dbl_phi    + (1-a) * tube_phi

    x,y,z = sphere_to_cart(theta, phi)

    new_data       = np.zeros(dbl_data.shape)
    new_data[:,:3] = dbl_data[:,:3]
    new_data[:,3]  = x
    new_data[:,4]  = y
    new_data[:,5]  = z

    return dim, new_data, num_vec, vec

def setup_ferro(v, cut_h=17):
    dim, data, num_vec, vec  = read_file("skyrm_tube.txt")
    dim, data = cut_data(dim, data, cut_h)

    for i in range(3):
        data[:,3+i] = v[i]
        
    return dim, data, num_vec, vec 

def insert_middle(filename, n_insert):
    dim, data, num_vec, vec  = read_file(filename)
    x   = data[:,0].reshape(dim, order="F")
    y   = data[:,1].reshape(dim, order="F")
    z   = data[:,2].reshape(dim, order="F")
    m_x = data[:,3].reshape(dim, order="F")
    m_y = data[:,4].reshape(dim, order="F")
    m_z = data[:,5].reshape(dim, order="F")

    dim[2] += n_insert
    print(dim)
    new_x = np.zeros(dim)
    new_y = np.zeros(dim)
    new_z = np.zeros(dim)
    new_m_x = np.zeros(dim)
    new_m_y = np.zeros(dim)
    new_m_z = np.zeros(dim)

    for i in range(dim[2]):
        new_x[:,:,i] = x[:,:,0]
        new_y[:,:,i] = y[:,:,0]
        new_z[:,:,i] = i

    middle = int(np.mean(z[0,0,:]))

    new_m_x[:,:,:middle]      = m_x[:,:,:middle]
    new_m_x[:,:,-(middle+2):] = m_x[:,:,-(middle+2):]
    new_m_y[:,:,:middle]      = m_y[:,:,:middle]
    new_m_y[:,:,-(middle+2):] = m_y[:,:,-(middle+2):]
    new_m_z[:,:,:middle]      = m_z[:,:,:middle]
    new_m_z[:,:,-(middle+2):] = m_z[:,:,-(middle+2):]

    for i in range(middle, dim[2]-(middle+2)):
        new_m_x[:,:,i] = m_x[:,:,middle]
        new_m_y[:,:,i] = m_y[:,:,middle]
        new_m_z[:,:,i] = m_z[:,:,middle]


    new_data = np.zeros((np.prod(dim), 6))
    new_data[:,0] = new_x.flatten(order="F")
    new_data[:,1] = new_y.flatten(order="F")
    new_data[:,2] = new_z.flatten(order="F")
    new_data[:,3] = new_m_x.flatten(order="F")
    new_data[:,4] = new_m_y.flatten(order="F")
    new_data[:,5] = new_m_z.flatten(order="F")

    return dim, new_data, num_vec, vec 

def top_skyrm(cut_h=-0.5):
    dim, data, num_vec, vec  = read_file("skyrm_tube.txt")
    data = ferro_below(data, 28.5)

    dim, data = cut_data(dim, data, cut_h)

    return dim, data, num_vec, vec
    
def center_inverted(cut_h=-0.5):
    dim, data, num_vec, vec  = read_file("skyrm_tube.txt")
    data = ferro_below(data, 30.0)

    dim, data = cut_data(dim, data, cut_h)

    pos = data[:,:3].copy()
    pos [:,0] -= np.mean(pos[:,0])
    pos [:,1] -= np.mean(pos[:,1])
    pos [:,2] -= np.mean(pos[:,2])

    pos = pos**2

    dist = np.sum(pos, axis=1)
    min_idx = np.argmin(dist)

    data[min_idx,3] = 0.0
    data[min_idx,4] = 0.0
    data[min_idx,5] = -1.0
    
    return dim, data, num_vec, vec


def write_file(fname, dim, data, num_vec, vec):
    print("saved {}".format(fname))
    with open(fname, mode="w") as f:
        f.write("# {} {} {}\n".format(dim[0], dim[1], dim[2]))

        for i in range(dim[0]*dim[1]*dim[2]):
            f.write("{:10} {:10} {:10} {:10} {:10} {:10}\n".format(
                data[i,0], data[i,1], data[i,2], data[i,3], data[i,4], data[i,5]))

        f.write("# {}\n".format(num_vec))
        for v in vec:
            f.write("{}".format(v))



dim, data, num_vec, vec  = read_file("bobber_large.txt")
data = ferro_below(data, 30)
dim, data = cut_data(dim, data, 8)
write_file("test.txt", dim, data, num_vec, vec)
