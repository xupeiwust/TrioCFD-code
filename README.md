
# TrioCFD

**TrioCFD** (previously named "Trio_U") is the Computational Fluid Dynamics (CFD) code
based on the TRUST platform.
There are different physical modules such as:
- Turbulence LES & RANS models,
- Front-Tracking,
- Radiation,
- ALE for fluid/structure interactions.

This software is OpenSource (BSD license).

Online documentation can be found at [https://triocfd-documentation.readthedocs.io/en/latest/](https://triocfd-documentation.readthedocs.io/en/latest/).

# **How to install TrioCFD-1.9.6 version ?**

### If TRUST-1.9.6 is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code#readme).

### Once TRUST installed, install TrioCFD-1.9.6 using one of these methods:

### **First method**
```bash
git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-1.9.6
cd TrioCFD-1.9.6
source PathToTRUST-1.9.6/env_TRUST.sh
baltik_build_configure -execute
make optim
make debug # TRUST should be built in debug mode
```

### **Second method**
```bash
wget ftp://ftp.cea.fr/pub/TRUST/TrioCFD/versions/v1.9.6/TrioCFD-1.9.6.tar.gz
tar xzf TrioCFD-1.9.6.tar.gz
mv TrioCFD TrioCFD-1.9.6
cd TrioCFD-1.9.6
source PathToTRUST-1.9.6/env_TRUST.sh
baltik_build_configure -execute
make optim
make debug # TRUST should be built in debug mode
```

# **How to install TrioCFD's development version ?**
**for developers and those interested in new features only.**

**Warning: "next" branch may not compile or some tests fail if important developments merged**

### If TRUST-next is not already installed, [please follow TRUST install instructions](https://github.com/cea-trust-platform/trust-code/tree/next#readme).
```bash
git clone https://github.com/cea-trust-platform/TrioCFD-code.git TrioCFD-next
cd TrioCFD-next
git checkout next
source PathToTRUST-next/env_TRUST.sh
baltik_build_configure -execute
make optim
make debug # TRUST should be built in debug mode
```
# **How to start ?**
```bash
source ./env_TrioCFD.sh
```

To check:
```bash
# All non-regression test cases in optimized mode:
make ctest_optim  # or make check_optim (slower)

# All non-regression test cases in debug mode:
make ctest_debug # or make check_debug (slower)

# A given non-regression test list (Some lists are available in ./share/testList/)
trust -ctest ./path/to/testList
```

To see documentation:
```bash
triocfd -index
```
