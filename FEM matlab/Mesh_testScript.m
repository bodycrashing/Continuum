
nc = 4;
no = 2;

%[nodes,elems] = Q4mesh_UnitCell(20,8,nc,no);
[nodes,elems,ndof] = MeshMaster('Q4',20,8,nc,no)

plot(nodes(:,1),nodes(:,2),'o')
hold on
test_x = reshape(nodes(:,1),nc+1,no+1);
test_y = reshape(nodes(:,2),nc+1,no+1);

plot(test_x,test_y)
hold on
plot(test_x',test_y')
axis equal
hold off
