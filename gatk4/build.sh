#!/bin/bash
BINARY_HOME=$PREFIX/bin
PACKAGE_HOME=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM
mkdir -p $PREFIX/bin
mkdir -p $PACKAGE_HOME

pip install --isolated biopython==1.70 bleach==1.5.0 cycler==0.10.0 enum34==1.1.6 h5py==2.7.1 html5lib==0.9999999 joblib==0.11 keras==2.2.0 markdown==2.6.9 matplotlib==2.1.0 pandas==0.21.0 patsy==0.4.1 protobuf==3.6.1 pymc3==3.1 pyparsing==2.2.0 pysam==0.13 python-dateutil==2.6.1 pytz==2017.3 pyvcf==0.6.8 pyyaml==3.12 scikit-learn==0.19.1 scipy==1.0.0 six==1.11.0 theano==1.0.4 tqdm==4.19.4 werkzeug==0.12.2

pip install gatkPythonPackageArchive.zip