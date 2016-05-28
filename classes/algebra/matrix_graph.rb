################################################################################
#
# Class MatrixGraph
#
# Author: Jan Rezac
# Date created: 2010-08-12
# License: Cuby license
# Description: Matrix-based representation of graphs, graph searching algorithms
# Status: Working, documented
#
################################################################################

module Algebra
class MatrixGraph < Matrix

	UNDEFINED = -1
	INFINITY = 1.0/0.0
	NOT_CONNECTED = INFINITY

	#=======================================================================
	# Constructor
	#=======================================================================

	def MatrixGraph.empty(size)
		return MatrixGraph.filled(size, size, NOT_CONNECTED)
	end

	#=======================================================================
	# Misc
	#=======================================================================

	def add_edge_undirected(i,j, weight)
		#: Add bidirectional edge, equivalent to:
		#| self[i,j] = weight
		#| self[j,i] = weight
		self[i,j] = weight
		self[j,i] = weight
	end

	def remove_edge_undirected(i,j)
		#: Remove undirected edge between i and j
		self[i,j] = NOT_CONNECTED
		self[j,i] = NOT_CONNECTED
	end

	def each_edge
		#: Iterates over all edges
		self.n.times{|i|
			self.n.times{|j|
				next if self[i,j] == NOT_CONNECTED
				yield(i,j)
			}
		}
		return nil
	end

	def each_edge_undirected
		#: Iterates over all edges of undirected graph 
		self.n.times{|i|
			i.times{|j|
				next if self[i,j] == NOT_CONNECTED
				yield(i,j)
			}
		}
		return nil
	end

	def edge_count_undirected # => integer
		#: Count edges in udirected graph
		count = 0
		self.n.times{|i|
			i.times{|j|
				count += 1 unless self[i,j] == NOT_CONNECTED
			}
		}
		return count
	end


	#=======================================================================
	# Best path searches
	#=======================================================================


	def dijkstra_best_path(source, target) # => Float, Array
		#: Return best distance and path from source to target. The path is array
		#: containing node numbers along the path from source to target.

		dist, previous = dijkstra(source, target)
		# Get path from the previous array
		path = get_path(source, target, previous)

		return [dist[target], path]
	end

	def all_possible_paths(source, target) #=> Array of Arrays
		#: Get all paths without loops connecting source and target and their lengts. The result is
		#: an array of [length, path] arrays. The list is ordered by path length.
	
		dist_s, previous_s = dijkstra(source)
		dist_t, previous_t = dijkstra(target)

		dists_and_paths = []

		# Each edge in graph
		each_edge{|i,j|
			# Assure no loops at the connection
			next if previous_s[i] == j
			next if previous_t[j] == i
		
			# Get distance & path
			dist = dist_s[i] + self[i,j] + dist_t[j]
			path = get_path(source, i, previous_s) + get_path_reverse(target, j, previous_t)

			# Skip paths with large loops
			next if path.size  != path.uniq.size

			save = true
			dists_and_paths.each{|d_a_p|
				if d_a_p[1] == path
					save = false
					break
				end
			}
			dists_and_paths << [dist, path] if save
		}

		dists_and_paths.sort!{|a, b| a[0] <=> b[0]}
		return dists_and_paths
	end

	def k_best_paths(source, target, k, one_for_each_dist = false)
		#: Get k best paths without loops connecting source and target and their lengts. The result is
		#: an array of [length, path] arrays. When one_for_each_dist parameter is set to true,
		#: only one (first found) path of each length is added to the list. The list is ordered by path length.

		if k == 1
			return [dijkstra_best_path(source, target)]
		end

		dist_s, previous_s = dijkstra(source)
		dist_t, previous_t = dijkstra(target)

		dists_and_paths = []

		# Each edge in graph
		each_edge{|i,j|
			# Assure no loops at the connection
			next if previous_s[i] == j
			next if previous_t[j] == i
		
			# Get distance
			dist = dist_s[i] + self[i,j] + dist_t[j]

			# Skip when it is worse than what we already have
			if dists_and_paths.size == k
				next if dists_and_paths.last[0] <= dist
			end

			# Skip if we already have path of the same length
			skip_this = false
			if one_for_each_dist
				dists_and_paths.each{|d_a_p|
					if dist == d_a_p[0]
						skip_this = true 
						break
					end
				}
			end
			next if skip_this

			# Get path
			path = get_path(source, i, previous_s) + get_path_reverse(target, j, previous_t)

			# Skip paths with large loops
			next if path.size  != path.uniq.size

			unique = true
			dists_and_paths.each{|d_a_p|
				if d_a_p[1] == path
					unique = false
					break
				end
			}

			if unique
				dists_and_paths << [dist, path] 
				dists_and_paths.sort!{|a, b| a[0] <=> b[0]}
				if dists_and_paths.size > k
					dists_and_paths.pop # Throw away last one
				end
			end
		}

		return dists_and_paths
	end

	#=======================================================================
	# Distance matrix algorithms
	#=======================================================================

	def floyd_warshall # => Matrix
		#: Plain Floyd-Warshall algorithm that return matrix of distances between all nodes
		n = self.n
		dists = Matrix.filled(n,n,INFINITY)
		n.times{|i| dists[i,i] = 0.0 } # Zero distance between self
		each_edge{|i,j| dists[i,j] = self[i,j]} # True distance between connected

		n.times{|k|
			n.times{|i|
				n.times{|j|
					tempdist = dists[i,k] + dists[k,j]
					if tempdist < dists[i,j]
						dists[i,j] = tempdist
					end
				}
			}
		}
		return dists
	end

	def floyd_warshall_saving_path #  => [Matrix, Array]
		#: Floyd-Warshall algorithm that return the distance matrix and data that can be used
		#: to reconstruct the paths (using -->floyd_warshall_path)
		n = self.n
		dists = Matrix.filled(n,n,INFINITY)
		n.times{|i| dists[i,i] = 0.0 } # Zero distance between self
		each_edge{|i,j| dists[i,j] = self[i,j]} # True distance between connected

		# Initialize array to save the path directions
		path_data = []
		n.times{|i| path_data[i] = [] }

		n.times{|k|
			n.times{|i|
				n.times{|j|
					tempdist = dists[i,k] + dists[k,j]
					if tempdist < dists[i,j]
						dists[i,j] = tempdist
						path_data[i][j] = k
					end
				}
			}
		}
		return [dists, path_data]
	end

	def floyd_warshall_path(i, j, path_data) # => Array
		#: Reconstructs shortest paths between nodes i and j using path data created by
		#: -->floyd_warshall_saving_path
		return [i] + floyd_warshall_path_priv(i, j, path_data) + [j]
	end

	def floyd_warshall_path_priv(i, j, path_data)
		#: Private method used for recursive generation of the paths
		intermediate = path_data[i][j]
		if intermediate == nil # No path, i and j connected by edge
			return []
		end
		return floyd_warshall_path_priv(i, intermediate, path_data) + [intermediate] + floyd_warshall_path_priv(intermediate, j, path_data)
	end
	private :floyd_warshall_path_priv

	#-----------------------------------------------------------------------

	def fw_with_equal_paths # => Matrix, Array
		#: Modification of the Floyd-Warshall algorithm that saves data that allow to recover all
		#: equivalent shortest paths
		n = self.n
		dists = Matrix.filled(n,n,INFINITY)
		n.times{|i| dists[i,i] = 0.0 } # Zero distance between self
		each_edge{|i,j| dists[i,j] = self[i,j]} # True distance between connected

		# Initialize array to save the path directions
		path_data = []
		n.times{|i| path_data[i] = [] }

		n.times{|k|
			n.times{|i|
				n.times{|j|
					tempdist = dists[i,k] + dists[k,j]
					if tempdist < dists[i,j]
						dists[i,j] = tempdist
						path_data[i][j] = [k]
					elsif tempdist == dists[i,j] && k != j && k != i
							if path_data[i][j]
								path_data[i][j] << k unless path_data[i][j].include?(k)
							end
					end
				}
			}
		}
		return [dists, path_data]
	end

	def fw_with_equal_paths_d_plus_1(prev_distmat) # => Matrix, Array
		#: Lists all distances and shortest paths that are longer by exactly 1.0 than previous result
		#: passed as distance matrix.
		n = self.n
		dists = Matrix.filled(n,n,INFINITY)

		# Initialize array to save the path directions
		path_data = []
		n.times{|i| path_data[i] = [] }

		n.times{|k|
			n.times{|i|
				n.times{|j|
					next if prev_distmat[i,j] == INFINITY
					tempdist = prev_distmat[i,k] + prev_distmat[k,j]
					if tempdist == prev_distmat[i,j] + 1.0
						dists[i,j] = tempdist
						if path_data[i][j] == nil
							path_data[i][j] = [k]
						else
							path_data[i][j] << k unless path_data[i][j].include?(k)
						end
					end
				}
			}
		}
		return [dists, path_data]
	end

	def fw_paths(i, j, path_data) # => Array
		#: Returns array of all equivalent paths between nodes i and j, using
		#: path data created by -->fw_with_equal_paths
		fullpaths = []
		fw_path_priv(i,j, path_data).each{|p|
			if p.size > 0 || i != j
				if p == []
					fullpaths << [i] + [j]
				else
					fullpaths << [i] + p + [j]
				end
			end
		}
		return fullpaths
	end

	def fw_paths_pdata_to_pmat(path_data) # => Array of Arrays of Arrays
		#: Builds 2D array (n x n) of arrays of equivalent paths between all nodes
		#: using path data created by -->fw_with_equal_paths
		n = self.n
		pmat = []
		n.times{|i|
			pmat[i] = []
			n.times{|j|
				pmat[i][j] = fw_paths(i, j, path_data)
			}
		}
		return pmat
	end

	def fw_paths_d_plus_1_pdata_to_pmat(path_data, pmat_prev)
		#: Builds 2D array (n x n) of arrays of equivalent (shortest + 1) paths between all nodes
		#: using path data created by -->fw_with_equal_paths_d_plus_1
		n = self.n
		pmat = []
		n.times{|i|
			pmat[i] = []
			n.times{|j|
				pmat[i][j] = fw_paths_d_plus_1_pmat(i, j, path_data, pmat_prev)
			}
		}
		return pmat
	end

	def fw_paths_d_plus_1_pmat(i, j, path_data, pmat_prev)
		retpaths = []
		return [] unless path_data[i][j]
		path_data[i][j].each{|k|
			pmat_prev[i][k].each{|ikpath|
				pmat_prev[k][j].each{|kjpath|
					np = ikpath[0..-2] + kjpath
					retpaths << np unless retpaths.include?(np)
				}
			}
		}
		return retpaths
	end

	def fw_paths_d_plus_1(i, j, path_data, path_data_prev)
		retpaths = []
		return [] unless path_data[i][j]
		path_data[i][j].each{|k|
			fw_paths(i, k, path_data_prev).each{|ikpath|
				fw_paths(k, j, path_data_prev).each{|kjpath|
					np = ikpath[0..-2] + kjpath
					retpaths << np unless retpaths.include?(np)
				}
			}
		}
		return retpaths
	end


	def fw_path_priv(i, j, path_data)
		intermediates = path_data[i][j]
		if intermediates == nil # No path, i and j connected by edge
			return [[]]
		end

		possibilities = []
		intermediates.each{|intermediate|
			fw_path_priv(i, intermediate, path_data).each{|seg_i|
				fw_path_priv(intermediate, j, path_data).each{|seg_j|
					tmp = seg_i + [intermediate] + seg_j
					possibilities << tmp unless possibilities.include?(tmp)
				}
			}
		}

		return possibilities
	end

	#=======================================================================
	# Ring search
	#=======================================================================

	def smallest_set_of_smallest_rings(size_limit = 9999999)

		# Algorithm from:
		# Lee, Kang, Cho and No, PNAS 106(41), 17355-17358, 2009

		### N_sssr is invalid in systems with disconnected rings
		n_sssr = edge_count_undirected - self.n + 1

		return [] if n_sssr == 0

		# Get path matrices
		distmat, path_data = fw_with_equal_paths
		distmat2, path_data2 = fw_with_equal_paths_d_plus_1(distmat)
		p = fw_paths_pdata_to_pmat(path_data)
		p2 = fw_paths_d_plus_1_pdata_to_pmat(path_data2, p)

		# Get all rings
		all_rings = []
		n.times{|i|
			i.times{|j|
				if p[i][j].size >= 2
					p[i][j].each_index{|a|
						a.times{|b|
							candidate = p[i][j][a] + p[i][j][b]
							orig_size = candidate.size
							candidate.uniq!
							next if candidate.size > size_limit
							if candidate.size == orig_size - 2
								# It's a ring!
								candidate.sort!
								all_rings << candidate unless all_rings.include?(candidate)
							end
						}
					}
				end
				if p[i][j].size == 1 && p2[i][j].size >= 1
					p2[i][j].each_index{|b|
						candidate = p[i][j][0] + p2[i][j][b]
						orig_size = candidate.size
						candidate.uniq!
						next if candidate.size > size_limit
						if candidate.size == orig_size - 2
							# It's a ring!
							candidate.sort!
							all_rings << candidate unless all_rings.include?(candidate)
						end
					}
				end
			}
		}

		# Sort all rings by size, lowest first
		all_rings.sort!{|a,b| a.size <=> b.size}

		# Get Smallest set of smallest rings (SSSR)
		sssr = []

		all_rings.each{|ring|
			valid = true
			sssr.each{|r_ok|
				if r_ok.size == (r_ok & ring).size
					valid = false
					break
				end
			}
			sssr << ring if valid
			#break if sssr.size == n_sssr
			### This does not work when there are multiple disconnected graphs in the system, n_sssr is then wrong
		}

		# Order elements in a ring to start at node with lowest index, continuing on direction of
		# its neighbor with lowest index
		sssr.each_index{|i|
			ring = sssr[i]
			new = []
			new << ring.shift
			while ring.size > 0
				ring.each_index{|j|
					if self[new.last, ring[j]] != NOT_CONNECTED
						new << ring[j]
						ring.delete_at(j)
						break
					end
				}
			end
			sssr[i] = new
		}

		return sssr

	end

	#=======================================================================
	# Maximum flow algorithm
	#=======================================================================
	
	def max_flow_edmonds_karp(source, sink) # => [Float, Matrix]
		#: Max flow algorithm of Edmonds and Karp.
		#: Return maximum possible flow through the graph and matrix
		#: of flows through all the edges. 

		# Adatpted from Java code at http://en.wikipedia.org/wiki/Edmonds%E2%80%93Karp_algorithm

		size = self.n
		flow = Matrix.zero(size,size)
		res_capacity = Matrix.zero(size,size)
		parent = Array.new(size)
		min_capacity = Array.new(size,0.0)

		max_flow = 0.0

		# Initialize residual capacity with maximum capacity
		each_index{|i,j|
			res_capacity[i,j] = self[i,j] if self[i,j] != INFINITY
		}

		while max_flow_bfs(source, sink, res_capacity, min_capacity, parent)
			max_flow += min_capacity[sink]
			v = sink
			while v != source
				u = parent[v]
				flow[u,v] += min_capacity[sink]
				flow[v,u] -= min_capacity[sink]
				res_capacity[u,v] -= min_capacity[sink]
				res_capacity[v,u] += min_capacity[sink]
				v = u
			end
		end

		return [max_flow, flow]
	end

	def max_flow_bfs(source, sink, res_capacity, min_capacity, parent)
		# Breadth-first search used by max_flow_edmonds_karp()
		
		size = self.n
		color = Array.new(size, :white)

		size.times{|i|
			min_capacity[i] = Float::MAX
		}

		queue = [source]
		color[source] = :gray

		while queue.size > 0
			v = queue.shift
			size.times{|u|
				if color[u] == :white && res_capacity[v,u] > 0
					min_capacity[u] = [min_capacity[v], res_capacity[v,u]].min
					parent[u] = v
					color[u] = :gray
					return true if u == sink
					queue << u
				end
			}
		end

		return false
	end
	private :max_flow_bfs
	
	#=======================================================================
	# Private methods - Dikstra an best path searches
	#=======================================================================
	
	def dijkstra(source, target = nil) #=> [Array, Array]
		#: Returns array of distances to visited nodes and array of parent in path for each node
		#: When target == nil, whole graph is traversed, otherwise the search i terminated wheb target is reached.

		# Algorithm based on http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

		n = self.n

		# Prepare temp. arrays
		in_q = Array.new(n, true)
		dist = Array.new(n, INFINITY)
		previous = Array.new(n, UNDEFINED)

		dist[source] = 0.0

		while in_q.include?(true)
			# get vertex in Q with smallest dist
			u = nil
			mindist = INFINITY
			n.times{|i|
				if in_q[i]
					if dist[i] < mindist
						mindist = dist[i]
						u = i
					end
				end
			}

			break if u == nil
			break if u == target if target != nil # otherwise whole graph is traversed

			in_q[u] = false	# remove u from Q

			# for each neighbor v of u  where v has not yet been removed from Q
			n.times{|v|
				next if ! in_q[v]
				next if self[u,v] == NOT_CONNECTED

				alt = dist[u] + self[u,v]
				if alt < dist[v]
					dist[v] = alt
					previous[v] = u
				end
			}
		end

		return dist, previous
	end
	private :dijkstra

	def get_path(source, target, previous) # => Array
		#: Convert raw results of dijkstra search to array listing nodes along the path from source to target.
		path = []
		u = target
		while previous[u] != UNDEFINED
			path.unshift u
			u = previous[u]
		end
		path.unshift source
		return path
	end
	private :get_path

	def get_path_reverse(source, target, previous) # => Array
		#: Equivalent of -->get_path, but with reverse ordering.
		path = []
		u = target
		while previous[u] != UNDEFINED
			path << u
			u = previous[u]
		end
		path << source
		return path
	end
	private :get_path_reverse
end
end

