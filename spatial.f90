sidelength = x
!current position for each body in all dimensions
!box variables have dimension for body
box_x(:) = r(1,:,0) / sidelength
box_y(:) = r(2,:,0) / sidelength
box_z(:) = r(3,:,0) / sidelength

if (box_x == i) then
    if (box_y == j) then
        if (box_z == k) then
            in_box(i,j,k) = true
        end if
    end if
end if


box_COM = 0
DO num = 1,9
    DO counter = 1, n
        !1st index is dimension,
        box_COM(:,num)) = box_COM(:,num) + r(:,counter,0)*M(counter)
    END DO
END DO




if (theta < 0.2) then !idk if 0.2 is a good guess

end if
