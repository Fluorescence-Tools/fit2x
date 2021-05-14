#ifndef FIT2X_CONDITIONALUPDATER_H
#define FIT2X_CONDITIONALUPDATER_H

#include <vector>
#include <memory>


class ConditionalUpdater{

private:

    /// false if model function needs to be updated
    bool _is_valid = false;

    std::vector<std::shared_ptr<ConditionalUpdater>> children = std::vector<std::shared_ptr<ConditionalUpdater>>();

public:

    bool get_is_valid() const {
        if(_is_valid){
            for(auto &v: children){
                if(!v->get_is_valid()){
                    return false;
                }
            }
        }
        return _is_valid;
    }

    void add_child(std::shared_ptr<ConditionalUpdater> child){
        children.emplace_back(child);
    }

protected:

    /// Used to set if the decay is 'valid'. A decay is invalid if the
    /// computed model, wres, etc. does not correspond to the input parameters.
    /// It should not be decided by outside if the decay is valid.
    void set_is_valid(bool v){
        _is_valid = v;
    }

    virtual void evaluate(){

    }

public:

    void update(){
        if(!get_is_valid()){
            for(auto &child:children){
                child->update();
            }
        }
        evaluate();
        set_is_valid(true);
    }

};

#endif //FIT2X_CONDITIONALUPDATER_H
