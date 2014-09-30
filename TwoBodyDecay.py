#!/usr/bin/env python
"""
Small program to calculate the kinematics of a two body decay of a mother
(particle 0) decaying into two daughters (particles 1 and 2). The program will
also perform a boost along the trajectory of particle 1, along the trajectory of
particle 2, and perpendicular to both particles.

One can either call this from the command line:
  ./TwoBodyDecay.py m0 p0 m1 m2

Or, one can import this into python and decay the particles interactively:
    decayParticle( m0, p0, m1, m2 )
"""

import math
import sys
import argparse

# ------------------------------------------------------------------------------
def getBeta(m, p):
    beta = math.sqrt( float(p)**2 / ( float(m)**2 + float(p)**2) )
    return beta

# ------------------------------------------------------------------------------
def getGamma(beta):
    return 1. / math.sqrt( 1 - float(beta)**2)

# ------------------------------------------------------------------------------
def getMassFromMomentum(p):
    return math.sqrt( p['e']**2 - p['px']**2 - p['py']**2 - p['pz']**2 )

# ------------------------------------------------------------------------------
def lorentzBoost(p0, beta, direction = 'x'):
    """
    Given an initial momentum 4-vector and a boost parameter beta, calculate the
    4-momentum of the particle in the boosted frame
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    gamma = getGamma(beta)
    if direction == 'x':
        return { 'e': gamma * ( p0['e']  - beta*p0['px'] )
               , 'px':gamma * ( p0['px'] - beta*p0['e']  )
               , 'py':p0['py']
               , 'pz':p0['pz']
               }
    if direction == 'y':
        return { 'e': gamma * ( p0['e']  - beta*p0['py'] )
               , 'px':p0['px']
               , 'py':gamma * ( p0['py'] - beta*p0['e']  )
               , 'pz':p0['pz']
               }
    if direction == 'z':
        return { 'e': gamma * ( p0['e']  - beta*p0['px'] )
               , 'px':p0['px']
               , 'py':p0['py']
               , 'pz':gamma * ( p0['pz'] - beta*p0['e']  )
               }
    return None

# ------------------------------------------------------------------------------
def restFrameMomenta(m0, m1, m2):
    """
    Calculate the momenta of the decay particles in the rest frame of the mother

    0 --> 1 2
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    p1 = 1/2. * math.sqrt( (m0**2 - (m1 + m2)**2))
    p2 = -p1

    return {'p1':p1, 'p2':p2}

# ------------------------------------------------------------------------------
def restFrameDecay(m0, m1, m2):
    daughter_momenta = restFrameMomenta(m0, m1, m2)
    p1 = { 'e': math.sqrt( m1**2 + daughter_momenta['p1']**2 )
         , 'px':daughter_momenta['p1']
         , 'py':0.
         , 'pz':0.
         }
    p2 = { 'e': math.sqrt( m2**2 + daughter_momenta['p2']**2 )
         , 'px':daughter_momenta['p2']
         , 'py':0.
         , 'pz':0.
         }

    return {'p1':p1, 'p2':p2}

# ------------------------------------------------------------------------------
def printDecay(p0, p1, p2):
    m0 = getMassFromMomentum(p0)
    m1 = getMassFromMomentum(p1)
    m2 = getMassFromMomentum(p2)

    p_0 = math.sqrt( p0['px']**2 + p0['py']**2 + p0['pz']**2 )
    p_1 = math.sqrt( p1['px']**2 + p1['py']**2 + p1['pz']**2 )
    p_2 = math.sqrt( p2['px']**2 + p2['py']**2 + p2['pz']**2 )

    print '---------------------------------------------'
    print '|    | %10s | %10s | %10s |' % ('mother', 'daughter 1', 'daughter 2')
    print '|-------------------------------------------|'
    print '| m  | %10.2f | %10.2f | %10.2f |' % ( m0, m1, m2 )
    print '| e  | %10.2f | %10.2f | %10.2f |' % ( p0['e'], p1['e'], p2['e'] )
    print '| p  | %10.2f | %10.2f | %10.2f |' % ( p_0, p_1, p_2 )
    print '| - - - - - - - - - - - - - - - - - - - - - |'
    print '| px | %10.2f | %10.2f | %10.2f |' % ( p0['px'], p1['px'], p2['px'] )
    print '| py | %10.2f | %10.2f | %10.2f |' % ( p0['py'], p1['py'], p2['py'] )
    print '| pz | %10.2f | %10.2f | %10.2f |' % ( p0['pz'], p1['pz'], p2['pz'] )
    print '---------------------------------------------'


# ------------------------------------------------------------------------------
def decayParticle(m0, p0, m1, m2):
    # first, get rest frame decay properties
    rest_frame_decay = restFrameDecay(m0, m1, m2)

    print ''
    print 'Rest frame decay'
    p0_rest = {'e':m0, 'px':0, 'py':0, 'pz':0}
    p1_rest = rest_frame_decay['p1']
    p2_rest = rest_frame_decay['p2']
    printDecay(p0_rest, p1_rest, p2_rest)

    if p0 == 0: return

    # calculate beta for boost calculations
    beta = getBeta(m0, p0)

    # Now, calculate boost in direction of particle 1
    print ''
    print 'Boost along particle 1 momentum'
    p0_boost_1 = lorentzBoost(p0_rest, -beta, 'x')
    p1_boost_1 = lorentzBoost(p1_rest, -beta, 'x')
    p2_boost_1 = lorentzBoost(p2_rest, -beta, 'x')
    printDecay(p0_boost_1, p1_boost_1, p2_boost_1)

    # Now, calculate boost in direction of particle 2
    print ''
    print 'Boost along particle 2 momentum'
    p0_boost_2 = lorentzBoost(p0_rest, beta, 'x')
    p1_boost_2 = lorentzBoost(p1_rest, beta, 'x')
    p2_boost_2 = lorentzBoost(p2_rest, beta, 'x')
    printDecay(p0_boost_2, p1_boost_2, p2_boost_2)

    # Now, calculate boost in direction perpendicular to  particle 1 and 2
    print ''
    print 'Boost perpendicular to particle 1 and particle 2'
    p0_boost_3 = lorentzBoost(p0_rest, -beta, 'y')
    p1_boost_3 = lorentzBoost(p1_rest, -beta, 'y')
    p2_boost_3 = lorentzBoost(p2_rest, -beta, 'y')
    printDecay(p0_boost_3, p1_boost_3, p2_boost_3)

# ==============================================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--m0', type=float, default=91)
    parser.add_argument('--p0', type=float, default=0)
    parser.add_argument('--m1', type=float, default=0.1)
    parser.add_argument('--m2', type=float, default=0.1)

    args = parser.parse_args()
    print args
    m0 = float(args.m0)
    p0 = float(args.p0)
    m1 = float(args.m1)
    m2 = float(args.m2)

    decayParticle( m0, p0, m1, m2 )
