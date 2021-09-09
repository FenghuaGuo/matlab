function  an1_2 = get_ang_from_dirs(dir1,dir2)

an1_2 = (180/pi)*real(acos(abs(sum(dir1.*dir2,1))));

end

function plot_dir_with_fof_glyphs(dir,val)

maxv = max(val);
mag = 5*maxv;

for ii = 1:length(dir)

    ps = mag*dir{1,ii};
    pe = -mag*dir{1,ii};
    plot3([ps(1);pe(1)],[ps(2);pe(2)],[ps(3);pe(3)])


end