# PKT-SMR-CPU-04 map lifetime

After daily reporting, the driver removes only flux, energy, temporary soil
redistribution, and zonal scratch maps. It retains all cross-day state maps:
Layer A/B and total storage, canopy storage, snow state, and saturation.
This bounds GRASS mapset growth without altering a model equation or step order.
