# ATLANTIS

R package for biomarker/depenedency relationship discovery (via learning random forests 
of conditional inference trees).

There are two packages contained in this directory which are both necessary
to run ATLANTIS:

    - partyMod: A fork from the party package which is used to actually
      build the model.  The fork was necessary to include tweaks to reduce 
      memory necessary for training with large numbers of input
      features.  (See https://cran.r-project.org/web/packages/party/index.html and
      "Torsten Hothorn, Peter Buehlmann, Sandrine Dudoit, Annette Molinaro and 
      Mark Van Der Laan (2006). Survival Ensembles. Biostatistics, 7(3),
      355--373.)

    - ATLANTIS: The wrapper around party which was used to learn the models
      for the paper "Defining a Cancer Dependency Map".  Includes some prefiltering
      of features and biased weights more heavily on sensitive lines.   See 
      the paper for more details on the specific settings used for each MDP
      class.
