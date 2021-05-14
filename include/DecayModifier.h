//
// Created by tpeulen on 4/2/2021.
//

#ifndef FIT2X_DECAYMODIFIER_H
#define FIT2X_DECAYMODIFIER_H

#include <memory>
#include "DecayCurve.h"

class DecayModifier{

public:

    virtual void resize(size_t n) = 0;

    virtual void add(DecayCurve& decay) = 0;

    virtual ~DecayModifier() = default;

};

#endif //FIT2X_DECAYMODIFIER_H
