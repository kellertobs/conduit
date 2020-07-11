function    dadr = ddr(a,r)

if length(r) == 1  % constant grid spacing, h = r
    dadr = diff(a,1,2)./r;
else
    dadr = diff(a,1,2)./diff(r,1,2);
end
