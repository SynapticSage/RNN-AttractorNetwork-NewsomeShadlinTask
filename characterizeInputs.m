function h = characterizeInputs(Itot,synapI,Iapp)
    
    h = figure(120);
    subaxis(3,2,1,1);
    imagesc(Itot'); title('I_{total}');
    subaxis(3,2,2,1);
    plot(mean(Itot,2)); title('\mu_{I_{total}}');
    subaxis(3,2,1,2);
    imagesc(Iapp'); title('I_{app}');
    subaxis(3,2,2,2);
    plot(mean(Iapp,2)); title('\mu_{I_{app}}');
    subaxis(3,2,1,3);
    imagesc(synapI'); title('I_{synap}');
    subaxis(3,2,2,3);
    plot(mean(synapI,2)); title('\mu_{I_{synap}}');

end