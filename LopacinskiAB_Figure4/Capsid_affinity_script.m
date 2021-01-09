%This script uses the Endres and Zlotnik, 2002 formalism to estimate median
%effective contact affinities from individual contact affinities of a
%picornaviral capsid.

close all hidden
clear all

sym = 5; %fivefold symmetry
Kcontact = 1e-3; %contact affinity for base parameter set

for num_contacts=1:5
    for num_subunits=2:12
        exponent=((num_subunits*num_contacts)/2)/(num_subunits-1);
        stat_term=(sym^(num_subunits-1)/num_subunits)^(-1/(num_subunits-1));
        Kapparent=stat_term*(1/Kcontact)^-exponent;
        if num_contacts <= num_subunits
            data(num_contacts,num_subunits-1)=Kapparent;
        else data(num_contacts,num_subunits-1)=nan;
        end
    end
end
heatmap(data)
median(median(data,'omitnan'),'omitnan')