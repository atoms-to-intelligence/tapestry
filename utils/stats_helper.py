# vim: tabstop=2 expandtab shiftwidth=2 softtabstop=8

# Input: list of numbers
#
# Outputs: list of (values, idx) in original, sorted by value
def get_sorted_indices(lst, invert=False):
  if not invert:
    sorted_lst = sorted(lst)
    sorted_idx = sorted(range(len(lst)), key = lambda k : lst[k])
  else:
    sorted_lst = sorted(lst, key = lambda k : -k)
    sorted_idx = sorted(range(len(lst)), key = lambda k : -lst[k])

  sorted_val_idx = ', '.join(
      [f'{sorted_lst[i]:.3f}:{sorted_idx[i]}' for i in range(len(lst))] )
  return sorted_val_idx


if __name__ == '__main__':
  lst = [ 30, 10, 40, 20, 70 ]
  print(get_sorted_indices(lst))
  print(get_sorted_indices(lst, invert=True))
