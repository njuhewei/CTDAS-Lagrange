standard_variables = { 'bio_flux_prior' : {'name'        : 'bio_flux_prior',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of carbon dioxide, terrestrial vegetation, not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_carbon_dioxide_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_opt' : {'name'          : 'bio_flux_opt',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of carbon dioxide, terrestrial biosphere , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_carbon_dioxide_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ocn_flux_prior' : {'name'        : 'ocn_flux_prior',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of carbon dioxide, open ocean , not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_carbon_dioxide_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ocn_flux_opt' : {'name'          : 'ocn_flux_opt',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of carbon dioxide, open ocean , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_carbon_dioxide_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fossil_flux_imp' : {'name'       : 'fossil_flux_imp',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of carbon dioxide, fossil fuel burning , imposed ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_carbon_dioxide_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'fire_flux_imp' : {'name'         : 'fire_flux_imp',\
                                         'units'         : 'mol m-2 s-1' ,\
                                         'long_name'     : 'Surface flux of carbon dioxide, biomass burning , imposed ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'surface_carbon_dioxide_mole_flux', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_prior_cov' : {'name'    : 'bio_flux_prior_cov',\
                                         'units'         : 'mol2 region-2 s-2' ,\
                                         'long_name'     : 'Covariance of surface flux of carbon dioxide, terrestrial vegetation , not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'bio_flux_opt_cov' : {'name'      : 'bio_flux_opt_cov',\
                                         'units'         : 'mol2 region-2 s-2' ,\
                                         'long_name'     : 'Covariance of surface flux of carbon dioxide, terrestrial vegetation , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ocn_flux_prior_cov' : {'name'    : 'ocn_flux_prior_cov',\
                                         'units'         : 'mol2 region-2 s-2' ,\
                                         'long_name'     : 'Covariance of surface flux of carbon dioxide, open ocean , not optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ocn_flux_opt_cov' : {'name'      : 'ocn_flux_opt_cov',\
                                         'units'         : 'mol2 region-2 s-2' ,\
                                         'long_name'     : 'Covariance of surface flux of carbon dioxide, open ocean , optimized ', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : '', \
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'decimal_date' :  {'name'         : 'decimal_date',\
                                         'units'         : 'years' ,\
                                         'long_name'     : 'dates and times', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'date', \
                                         'dims'          : (), \
                                         'dtype'         : 'double', \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'date' :         {'name'          : 'date',\
                                         'units'         : 'days since 2000-01-01 00:00:00 UTC' ,\
                                         'long_name'     : 'UTC dates and times', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'standard_name' : 'date', \
                                         'dims'          : (), \
                                         'dtype'         : 'double', \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'idate' :        {'name'          : 'idate',\
                                         'units'         : 'yyyy MM dd hh mm ss ' ,\
                                         'long_name'     : 'integer components of date and time', \
                                         'standard_name' : 'calendar_components', \
                                         'comment'       : 'time-interval average, centered on times in the date axis', \
                                         'dims'          : (), \
                                         'dtype'         : 'int', \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'latitude' :     {'name'          : 'latitude',\
                                         'units'         : 'degrees_north ' ,\
                                         'long_name'     : 'latitude', \
                                         'standard_name' : 'latitude', \
                                         'comment'       : 'center of interval',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'longitude' :     {'name'         : 'longitude',\
                                         'units'         : 'degrees_east ' ,\
                                         'long_name'     : 'longitude', \
                                         'standard_name' : 'longitude', \
                                         'comment'       : 'center of interval',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'height' :        {'name'         : 'height',\
                                         'units'         : 'masl ' ,\
                                         'long_name'     : 'height_above_ground_level', \
                                         'standard_name' : 'height_above_ground_level', \
                                         'comment'       : 'value is meters above sea level',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'co2' :           {'name'         : 'co2',\
                                         'units'         : 'micromol mol-1 ' ,\
                                         'long_name'     : 'mole_fraction_of_carbon_dioxide_in_air', \
                                         'standard_name' : 'mole_fraction_of_carbon_dioxide_in_air', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'meanstate' :     {'name'         : 'statevectormean',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'mean_value_of_state_vector', \
                                         'standard_name' : 'mean_value_of_state_vector', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'ensemblestate':  {'name'         : 'statevectorensemble',\
                                         'units'         : 'unitless' ,\
                                         'long_name'     : 'ensemble_value_of_state_vector', \
                                         'standard_name' : 'ensemble_value_of_state_vector', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                       'unknown' :      {'name'          : '',\
                                         'units'         : '' ,\
                                         'long_name'     : '', \
                                         'standard_name' : '', \
                                         'comment'       : '',\
                                         'dims'          : (), \
                                         'values'        : [], \
                                         'count'         : 0 \
                                        } , \
                     }




