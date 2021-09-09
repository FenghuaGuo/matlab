figure; 

for ii= 1:9
    subplot(3,3,ii)
    plot_dir_with_fof_glyphs(dir_ini{1,ii},val_ini{1,ii},fof(:,ii),nd)
end

figure; 

for ii= 1:9
    subplot(3,3,ii)
    plot_dir_with_fof_glyphs(dir_ini{1,end-ii+1},val_ini{1,end-ii+1},fof(:,end-ii+1),nd)
end