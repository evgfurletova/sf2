#norootforbuild

Name:		sf2
Summary:	SufPref program to calculate P-value
Version:	2.0.0
Release:	1%{?dist}
License:	MIT
BuildRequires:	cmake >= 2.6
BuildRequires:	python-devel >= 2.6
BuildRequires:	gcc
BuildRequires:	gcc-c++
Source:		%name-%version.tar.gz


%description
A novel algorithm SufPref that computes an exact P-value for Hidden 
Markov models (HMM). The algorithm inductively traverses a specific 
data structure - overlap graph. The nodes of the graph are associated 
with the overlaps of words from H. The edges are associated to the 
prefix and suffix relations between overlaps. An originality of our 
data structure is that pattern H need not be explicitly represented 
in nodes or leaves. The algorithm relies on the Cartesian product of 
the overlap graph and the graph of HMM states. The gain in size of 
SufPref data structure leads to significant improvements in space and 
time complexity compared to existent algorithms.


%prep
%setup -q 


%build
cmake -DCMAKE_INSTALL_PREFIX=$RPM_BUILD_ROOT/%_prefix -DCMAKE_BUILD_TYPE=Release
make


%install
make install

%clean
rm -rf $RPM_BUILD_ROOT

%package sufpref
Summary:	SufPref executable

%description sufpref
The command line implementation of SufPref algorithm.

A novel algorithm SufPref that computes an exact P-value for Hidden 
Markov models (HMM). The algorithm inductively traverses a specific 
data structure - overlap graph. The nodes of the graph are associated 
with the overlaps of words from H. The edges are associated to the 
prefix and suffix relations between overlaps. An originality of our 
data structure is that pattern H need not be explicitly represented 
in nodes or leaves. The algorithm relies on the Cartesian product of 
the overlap graph and the graph of HMM states. The gain in size of 
SufPref data structure leads to significant improvements in space and 
time complexity compared to existent algorithms.

%files sufpref
%defattr(-,root,root)
%_bindir/sufpref

%package -n python-sf
Summary:	SufPref python API module

%description -n python-sf
The Python module implementing SufPref algorithm.

A novel algorithm SufPref that computes an exact P-value for Hidden 
Markov models (HMM). The algorithm inductively traverses a specific 
data structure - overlap graph. The nodes of the graph are associated 
with the overlaps of words from H. The edges are associated to the 
prefix and suffix relations between overlaps. An originality of our 
data structure is that pattern H need not be explicitly represented 
in nodes or leaves. The algorithm relies on the Cartesian product of 
the overlap graph and the graph of HMM states. The gain in size of 
SufPref data structure leads to significant improvements in space and 
time complexity compared to existent algorithms.

%files -n python-sf
%defattr(-,root,root)
%python_sitearch/sf.so

%package -n python-sf-doc
BuildArch:	noarch
Summary:	Documentation on SufPref Python 

%description -n python-sf-doc
SufPref Python API developer reference

%files -n python-sf-doc
%defattr(-,root,root)
%_docdir/%name/py_module.html
