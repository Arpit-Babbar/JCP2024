using Suppressor
using Trixi2Vtk
using ReadVTK
last_file_index = 29900
step = 100
n_files = Int(last_file_index / step)
n_elements_list = zeros(n_files)
for i in 1:n_files
unpadded_index = Int(i*step)
if unpadded_index < 1000000
	total_digits = 6
else
	total_digits = 7
end
file_index = lpad(unpadded_index, total_digits, "0")
local filename = "solution_$file_index"
if !(isfile("$filename.vtu"))
    @suppress trixi2vtk("$filename.h5")
end
n_elements = VTKFile("$filename.vtu").n_cells
println(n_elements)
n_elements_list[i] = n_elements
end
