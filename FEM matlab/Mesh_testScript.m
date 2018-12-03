
nc = 2;
no = 1;

%[nodes,elems] = Q4mesh_UnitCell(20,8,nc,no);
[nodes,elems,ndof] = MeshMaster('LSTiso',1,0.4,nc,no);

plot(nodes(:,1),nodes(:,2),'o')
hold on
%test_x = reshape(nodes(:,1),nc+1,no+1);
%test_y = reshape(nodes(:,2),nc+1,no+1);

elems = elems(:,[1,4,2,5,3,6]);

for i = 1:length(elems)

    x = nodes(elems(i,:),1);
    y = nodes(elems(i,:),2);
    
    x = [x;x(1)];
    y = [y;y(1)];
    
plot(x,y,'k')
plot(x',y','k')



end
axis equal
hold off