function plot_dir_with_fof_glyphs(dir,val,fod,dirs)

plot_odf_Chantal([fod;fod],[dirs;-dirs])
hold on
maxv = max(val);
mag = 1.5*maxv;

for ii = 1:length(dir)
%     mag = 3*val(ii);
    ps = mag*dir(:,ii);
    pe = -mag*dir(:,ii);
    plot3([ps(1);pe(1)],[ps(2);pe(2)],[ps(3);pe(3)])
end

hold off


end

