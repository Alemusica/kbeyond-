#include "c74_max.h"
#include <algorithm>
#include <memory>

using namespace c74::max;

namespace {

constexpr double kMinPhiweight = 0.0;
constexpr double kMaxPhiweight = 1.0;

t_class* s_kbeyond_class = nullptr;

struct kbeyond_engine {
    double phiweight {0.25};

    void update_householder() {}
};

struct t_kbeyond {
    t_pxobject ob;
    double phiweight {0.25};
    std::unique_ptr<kbeyond_engine> _kbeyond;
};

void* kbeyond_new(t_symbol*, long, t_atom*);
void kbeyond_free(t_kbeyond*);
void kbeyond_assist(t_kbeyond*, void*, long, long, char*);
void kbeyond_attr_set_phiweight(t_kbeyond* x, void* attr, long argc, t_atom* argv);

void* kbeyond_new(t_symbol*, long, t_atom*) {
    auto* x = reinterpret_cast<t_kbeyond*>(object_alloc(s_kbeyond_class));
    if (!x)
        return nullptr;

    dsp_setup(reinterpret_cast<t_pxobject*>(x), 1);
    x->_kbeyond = std::make_unique<kbeyond_engine>();

    atom args[1];
    atom_setfloat(args, x->phiweight);
    kbeyond_attr_set_phiweight(x, nullptr, 1, args);

    return x;
}

void kbeyond_free(t_kbeyond* x) {
    dsp_free(reinterpret_cast<t_pxobject*>(x));
}

void kbeyond_assist(t_kbeyond*, void*, long, long, char*) {}

void kbeyond_attr_set_phiweight(t_kbeyond* x, void*, long argc, t_atom* argv) {
    if (!x)
        return;

    double value = x->phiweight;
    if (argc > 0 && argv) {
        if (atom_gettype(argv) == A_LONG)
            value = static_cast<double>(atom_getlong(argv));
        else if (atom_gettype(argv) == A_FLOAT)
            value = atom_getfloat(argv);
    }

    value = std::clamp(value, kMinPhiweight, kMaxPhiweight);

    x->phiweight = value;
    if (x->_kbeyond) {
        x->_kbeyond->phiweight = value;
        x->_kbeyond->update_householder();
    }
}

} // namespace

int C74_EXPORT main() {
    t_class* c = class_new("kbeyond~", reinterpret_cast<method>(kbeyond_new), reinterpret_cast<method>(kbeyond_free),
        sizeof(t_kbeyond), reinterpret_cast<method>(nullptr), A_GIMME, 0);

    class_addmethod(c, reinterpret_cast<method>(kbeyond_assist), "assist", A_CANT, 0);

    CLASS_ATTR_DOUBLE(c, "phiweight", 0, t_kbeyond, phiweight);
    CLASS_ATTR_ACCESSORS(c, "phiweight", NULL, kbeyond_attr_set_phiweight);
    CLASS_ATTR_FILTER_CLIP(c, "phiweight", kMinPhiweight, kMaxPhiweight);

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    s_kbeyond_class = c;

    return 0;
}
