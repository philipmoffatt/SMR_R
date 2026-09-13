# PKT-SMR-CPU-03 single-watershed reporting

When `wshed_res_properties.ini` has one nonblank record, the GRASS watershed
MASK makes `r.univar -g map=<name>` a direct sum over the same cells. The
reporting loop uses that sum instead of creating a temporary zonal raster and
then running `r.stats`. More than one record retains the legacy zonal route.

This changes reporting processes only; no state map expression changes.
