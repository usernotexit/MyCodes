[x, t] = readObj('../obj/Balls.obj');

figure; subplot(131); trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis equal; axis off; title('coarse mesh');
[x2,t2] = loop(x,t);
%[v,e,t2,T] = loop2(x,t); x2 = [v;e];
subplot(132); trimesh(t2, x2(:,1), x2(:,2), x2(:,3), 'edgecolor', 'k'); axis equal; axis off; title('fine mesh');
[x3,t3] = loop(x2,t2);
subplot(133); trimesh(t3, x3(:,1), x3(:,2), x3(:,3), 'edgecolor', 'k'); axis equal; axis off; title('finer mesh');           