Convolutions
============

make figure of different convolution modes describe what to use when.

Pile-up correction
==================
Typical time-correlated single photon counting electronics detects can only
detect a single photon per excitation pulse. Thus, for a certain amount of time
after detecting a photon the detection system is unable to register photons.
Consequently, at higher count rates not all photons are registered. The count
rate depends on the time since excitation and decays (in a fluorescence experiment)
typically exponentially. Hence, the most photons are detected right after exciting
the sample. Subsequent photons are detected less likely. This results in an
apparent "pile-up" of photons in the beginning of the decay curve. This pile-up
is in fact the loss of photons that would have been detected at later times.

.. plot:: plots/pile_up.py

    Effect of pile-up on the shape of fluorescence decays

The effect of pile-up is well known and there are models that can compute the
effect of pile-up on experimental decay curves :cite:´coates_correction_1968´.
To compute the effect of pile-up the repetition rate of the laser, the dead time
of the detection system, and the total measurement time needs to be known.
Traditionally, the data can be corrected for pile-up :cite:´coates_correction_1968´.
However, for analysis, the model function can be perturbed to mirror the effect
of pile-up on the experimental data and to preserve the counting statistics.


Differential non linearities
============================

