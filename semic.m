% Eq.(3.24) Here c is half of that in Eq.(3.24).
function dy = semic(~,y0)
    global ntra q c flag
    
    y = reshape(y0,3,ntra);
    r1 = abs(y(1,:)).^2;
    r2 = abs(y(2,:)).^2;
    r3 = abs(y(3,:)).^2;
    dy = zeros(3,ntra);    % a column vector
    
    dy(1,:) = q*y(1,:)+c*((r1+r2-r3).*y(1,:)+flag*y(2,:).^2.*conj(y(3,:)));
    dy(2,:) = c*((r1+r3).*y(2,:)+flag*2*y(1,:).*y(3,:).*conj(y(2,:)));
    dy(3,:) = q*y(3,:)+c*((r3+r2-r1).*y(3,:)+flag*y(2,:).^2.*conj(y(1,:)));
    
    dy = -1i*dy;
    dy = dy(:);
end

