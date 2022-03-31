#!/usr/bin/env python

from taurunner.cross_sections import CrossSections
from taurunner.body import Body
from taurunner.body.earth import construct_earth
from taurunner.utils.make_tracks import make_tracks
from taurunner.utils import make_propagator
import taurunner
import pyrex
from pyrex.internal_functions import get_from_enum
from pyrex.particle import Interaction, GQRSInteraction, CTWInteraction

import numpy as np
import sys

def calc_distE(eini, th_nadir, phi, rand, xs, tau_prop, Earth, ptype):
    particle = taurunner.particle.Particle(ID=ptype, energy=eini, position=0., rand=rand, xs=xs,
                                       proposal_propagator=tau_prop, secondaries=True, no_losses=False) 
    thetas_nadir = np.ones(1)*th_nadir
    tracks = make_tracks(thetas_nadir)
    mytrack = tracks[thetas_nadir[0]]

    E0 = (particle.energy/taurunner.utils.units.TeV)*1e-3
    particle.Interact('CC', Earth, mytrack, proton_fraction=Earth.get_proton_fraction(mytrack.x_to_r(particle.position)))
    E1 = (particle.energy/taurunner.utils.units.TeV)*1e-3
    y = 1.-(E1/E0)

    sec = particle.PropagateChargedLepton(Earth, mytrack)
    particles = sec.particles
    decay_products = [p for i,p in zip(range(max(len(particles)-3,0),len(particles)), particles[-3:]) if int(p.type) <= 1000000001]
    if len(decay_products)>0:   
        dist = particle.chargedposition*1e3
        decayE = (particle.energy/taurunner.utils.units.TeV)*1e-3
        return dist, decayE, y
    else:
        return None, None, None

primary_E = 1e18
outer_radius = 1e5
outer_z = 1e5
inner_radius = 15000.
inner_z = 2800.
seednum = int(sys.argv[1])
primary_E_str = str(sys.argv[2])
primary_E = float(sys.argv[2])

xs_model = 'CSMS'

Earth    = construct_earth(layers=[(4., 1.0)])
xs       = CrossSections(xs_model)
rand     = np.random.RandomState()
tau_prop = make_propagator(16, Earth)

earth_pyrex = pyrex.earth_model.PREM()

generator = pyrex.CylindricalGenerator(
    dr=outer_radius, dz=outer_z,
    energy=primary_E*1e-9,
    shadow=False,
    earth_model=earth_pyrex,
    source="cosmogenic",
)

vertexes = []
dirs = []
Es = []
ptypes = []
gens = []
ys = []
sweights = []
iweights = []
itypes = []
for i in range(0,10000000, 10):
    
    primary_evt = generator.create_event()
    primary_p = primary_evt.roots[0]

    primary_vertex = primary_p.vertex
    primary_dir = primary_p.direction
    primary_type = primary_p.id
    primary_weights = generator.get_weights(primary_p)
    p_sweight = primary_weights[0]
    p_iweight = primary_weights[1]
    primary_y = primary_p.interaction.choose_inelasticity()

    primary_kind = primary_p.interaction.kind

    dist_to_center_primary = np.sqrt(primary_vertex[0]**2+primary_vertex[1]**2)

    if dist_to_center_primary < inner_radius and primary_vertex[2]>inner_z*-1:
        print("primary in sim volume")
        print(primary_vertex, primary_dir, primary_E)
        vertexes.append(primary_vertex)
        dirs.append(primary_dir)
        Es.append(primary_E*1e-15)
        ptypes.append(primary_type.value)
        gens.append(i)
        ys.append(y)
        sweights.append(p_sweight)
        iweights.append(p_iweight)
        if primary_kind == Interaction.Type.charged_current:
            itypes.append('CC')
        if primary_kind == Interaction.Type.neutral_current:
            itypes.append('NC')

    if (primary_type == pyrex.Particle.Type.tau_antineutrino or primary_type == pyrex.Particle.Type.tau_neutrino) and primary_kind == Interaction.Type.charged_current:
        norm_factor = np.sqrt(primary_dir[0]**2+primary_dir[1]**2+primary_dir[2]**2)
        normed_dir = primary_dir/norm_factor

        th = np.arccos(normed_dir[2])
#        th = 0.0#This is temporary!
        phi = np.arctan(normed_dir[1]/normed_dir[0])
        th_nadir = np.pi-th
        
        th_nadir = np.degrees(th_nadir)
        phi = np.degrees(phi)
        try:
            test_L, test_Edecay, y = calc_distE(primary_E, th_nadir, phi, rand, xs, tau_prop, Earth, primary_type.value)

            if test_L != None:
                secondary_vertex = primary_vertex+(test_L*primary_dir)
                dist_to_center = np.sqrt(secondary_vertex[0]**2+secondary_vertex[1]**2)
                
                dist_to_center_primary = np.sqrt(primary_vertex[0]**2+primary_vertex[1]**2)
                
                if dist_to_center < inner_radius and secondary_vertex[2]>inner_z*-1 and secondary_vertex[2]<0.:
                    print("secondary in sim volume")
                    print(secondary_vertex, normed_dir, test_L, test_Edecay, dist_to_center, inner_radius)
                    secondary_p = pyrex.Particle(particle_id=pyrex.Particle.Type.tau_antineutrino,
                                                 vertex = primary_vertex,
                                                 direction = primary_dir,
                                                 energy=test_Edecay,
                                                 interaction_model=CTWInteraction,
                                                 #interaction_type=Interaction.Type.charged_current)
                                                 interaction_type=None)
                    ccnc = secondary_p.interaction.kind
                    print("type of interaction is", ccnc)
                    y2 = secondary_p.interaction.choose_inelasticity()
                    print(type(primary_p), type(secondary_p))
                    s_iweight = generator.get_weights(secondary_p)[1]
                    vertexes.append(secondary_vertex)
                    dirs.append(normed_dir)
                    Es.append(test_Edecay)
                    gens.append(i+1)
                    #ys.append(None)
                    ys.append(y2)
                    sweights.append(p_sweight)
                    iweights.append(s_iweight)
                    itypes.append('CC')
                    if primary_type == pyrex.Particle.Type.tau_antineutrino:
                        ptypes.append(pyrex.Particle.Type.antitau.value)
                    elif primary_type == pyrex.Particle.Type.tau_neutrino:
                        ptypes.append(pyrex.Particle.Type.tau.value)
                else:
                    pass
            else:
                pass

        except IndexError:
            print("Tau did not decay in simulated volume")
            pass

    else:
        pass
#    else:
#        if dist_to_center_primary < inner_radius and primary_vertex[2]>inner_z*-1:
#            print("primary in sim volume")
#            print(primary_vertex, primary_dir, primary_E)
#            ccnc = primary_p.interaction.kind
#            vertexes.append(primary_vertex)
#            dirs.append(primary_dir)
#            Es.append(primary_E*1e-15)
#            ptypes.append(primary_type.value)
#            gens.append(i)
#            ys.append(primary_y)
#            sweights.append(p_sweight)
#            iweights.append(p_iweight)
#
#        else:
#            pass

vertexes = np.array(vertexes)
dirs = np.array(dirs)
Es = np.array(Es)
ys = np.array(ys)
sweights = np.array(sweights)
iweights = np.array(iweights)
itypes = np.array(itypes)

vert_x = [v[0] for v in vertexes]
vert_y = [v[1] for v in vertexes]
vert_z = [v[2] for v in vertexes]

dir_x = [d[0] for d in dirs]
dir_y = [d[1] for d in dirs]
dir_z = [d[2] for d in dirs]

outarr = np.recarray((len(vertexes),), dtype=[('vert_x', float), ('vert_y', float), ('vert_z', float), ('dir_x', float), ('dir_y', float), ('dir_z', float), ('energy', float), ('y', float) ,('ptype', int),('sweight', float), ('iweight', float), ('itype', str), ('generation', int)])

outarr['vert_x'] = vert_x
outarr['vert_y'] = vert_y
outarr['vert_z'] = vert_z
outarr['dir_x'] = dir_x
outarr['dir_y'] = dir_y
outarr['dir_z'] = dir_z
outarr['energy'] = Es
outarr['y'] =  ys
outarr['ptype'] = ptypes
outarr['generation'] = gens
outarr['sweight'] = sweights
outarr['iweight'] = iweights
outarr['itype'] = itypes

print(outarr)
np.save('/data/user/wluszczak/ara/secondaries/%s/dummy_evts_%s_%s.npy'%(primary_E_str, seednum, primary_E_str), outarr)
