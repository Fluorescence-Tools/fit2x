Glossary
========


.. glossary::

    CLSM
        confocal laser scanning microscopy

    MFD (Multiparameter Fluorescence Detection)
        A MFD experiments is a time-resolved fluorescence experiment which probes
        the absorption and fluorescence, the fluorescence quantum yield, the
        fluorescence lifetime, and the anisotropy of the studied chromophores
        simultaneously (see :cite:`Kuehnemuth2001`)

    IRF
        IRF stands for instrument response function. In time-resolved fluorescence
        measurements the IRF is the temporal response of the fluorescence spectrometer
        to a delta-pulse. Suppose a initially sharp pulse defines the time of
        excitation / triggers the laser, then recorded response of the fluorescence
        spectrometer is broadened due to: (1) the temporal response of the exciting
        light source, (2) the temporal dispersion due to the optics of the instrument,
        (3) the delay of the light within the sample, and (4) the response of the
        detector. As the most intuitive contribution to the IRF is the excitation
        profile, the IRF is sometimes called 'lamp function'. The IRF is typically
        recorded by minimising the contribution of (3), e.g., by measuring the response
        of the instrument using a scattering sample, or a short lived dye.

    Time-tagged time resolved (TTTR)
        TTTR stands for time tagged time-resolved data or experiments. In
        TTTR-datasets the events, e.g., the detection of a photon, are tagged by
        a detection channel number. Moreover, the recording clock usually registers the
        events with a high time resolution of a few picoseconds. For long recording
        times of the detected events, a coarse and a fine clock are combined. The
        fine clock measures the time of the events relative to the coarse
        clock with a high time resolution. The time of the coarse and the fine
        clock is usually called macro and micro time, respectively.

    Time correlated single photon counting (TCSPC)
        Time correlated single photon counting (TCPSC) is a technique to measure
        light intensities with picosecond resolution. Its main application is the
        detection of fluorescent light. A pulsed light source excites a fluorescent
        sample. A single photon detector records the emitted fluorescence photons.
        Thus, per excitation cycle, only a single photon is detected. Fast detection
        electronics records the time between the excitation pulse and the detection
        of the fluorescence photon. A histogram accumulates multiple detected photons
        to yield a time-resolved fluorescence intensity decay.

    SWIG
        SWIG is a software development tool that connects programs written in C
        and C++ with a variety of high-level programming languages. SWIG can be
        used with different types of target languages including common scripting
        languages such as Javascript, Perl, PHP, Python, Tcl and Ruby and
        non-scripting languages such as C#, D, Go language, Java, Octave, and R.
        SWIG is free software and the code that SWIG generates is compatible with
        both commercial and non-commercial projects. ``fit2x`` is C/C++ based
        to provide the capability for a broad variety of languages to interface
        its provided functionality.

    Scatter fraction
        The scatter fraction :math:`gamma` is defined by the number of photons
        that

    Anisotropy
        The steady-state anisotropy :math:`r_G` in the detection channel :math:`G`
        is formally given by the fluorescence intensity weighted integral of the
        time-resolved anisotropy.

        :math:`r_G=\int F_G(t) \cdot r(t) dt \cdot \frac{1}{\int F_G(t) dt}`

        where the time-resolved anisotropy is defined by unperturbed the fluorescence
        intensities of an ideal detection system.

        :math:`r_G(t)=\frac{F_{G,p}(t)-F_{G,s}(t)}{F_{G,p}(t)+2F_{G,s}(t)}`

        Through out ``fit2x`` two distinct anisotropies are computed: (1)
        background corrected anisotropies, and (2) anisotropies not accounting for
        the background. In single-molecule experiments the background is mainly
        scattered light (Raman scattering). The uncorrected anisotropy (without
        background correction) is computed by:

        :math:`r = (S_p - g \cdot S_s) / (S_p \cdot (1 - 3 \cdot l_2) + (2 - 3 \cdot l_1) \cdot g \cdot Ss)`

        where :math:`S_p` is the signal in the parallel (German: parallel=p) detection
        channel, :math`S_s` the signal in the perpendicular decection channel
        (German: senkrecht=s), :math:`g` is the g-factor, :math:`l_1` and
        :math:`l_2` are factor mixing that determine the mixing of the parallel
        and perpendicular detection channel, respectively :cite:`koshioka_time-dependent_1995`.

        The scatter corrected steady-state anisotropy is computed using the scatter /
        background corrected signals parallel :math:`F_p = (S_p - \gamma \cdot B_p) / (1. - \gamma)`
        and perpendicular :math:`F_s = (S_s - \gamma \cdot B_s) / (1. - \gamma)`
        fluorescence intensity.
        :math:`r = (F_p - g \cdot F_s) / (F_p \cdot (1 - 3 \cdot l_2) + (2 - 3 \cdot l_1) \cdot g \cdot F_s)`
        The scatter corrected and anisotropy not corrected for scatter are computed
        by most fits of ``fit2x``.

    Jordi-format
        In the Jordi format is a format for fluorescence decays. In the Jordi
        format fluorescence decays are stacked in a one dimensional array.
        In a typical polarization resolved Jordi file the first decay is
        the parallel and the subsequent decay is the perpendicular decay. In the
        Jordi format both decays must have the same length, i.e., the same number
        of micro time counting channels.