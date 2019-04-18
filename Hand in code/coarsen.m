function Rc=coarsen(R,m)
    mc = (m-1)/2;
    RR = zeros(m,m);
    for ii = 1:m
        RR(m-(ii-1),:) = R(1+(ii-1)*m:ii*m);
    end
    RR = RR(2:m-1,2:m-1);
%     size(RR)
    Rc = zeros(mc^2,1)+10;
    %hc = 1/(mc+1);
    for ii = 1:mc
        for jj = 1:mc
            Rc(jj+(ii-1)*mc) = RR((m-2)-2*(ii-1),1+2*(jj-1));
        end
    end
end



