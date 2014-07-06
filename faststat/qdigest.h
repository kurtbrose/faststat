/*
An implementation of the QDigest algorithm for floating point numbers.

(The internal data structure is kept in terms of the bits of the floating point.)
*/

#ifdef _WIN32
#define inline __inline
#endif

#include <Python.h>

typedef struct {
    short next;
    short padddd;
    unsigned int min;
    unsigned long long count;
} Qdigest_node;


typedef struct Qdigest {
    unsigned long long n;
    short k;
    short free_head;
    short free_tail;
    short num_free;
    short new_head;
    short new_tail;
    short generations[32];
    Qdigest_node *nodes; // TODO: flexible array?  (does MSVC support this feature?)
} Qdigest;


Qdigest* qdigest_new(short size) {
    Qdigest *q;
    Qdigest_node *nodes;
    int i;
    // allocate memory
    q = PyMem_Malloc(sizeof(Qdigest));
    nodes = PyMem_Malloc(sizeof(Qdigest_node) * size);
    q->k = size / 4;
    q->n = 0;
    q->nodes = nodes;
    // all nodes are free to start with
    q->free_head = 0;
    q->free_tail = size - 1;
    for(i=0; i < size - 1; i++) {
        nodes[i].next = i + 1;
    }
    nodes[size - 1].next = -1;
    q->num_free = size;
    // all generations are empty
    for(i=0; i < 32; i++) {
        q->generations[i] = -1;
    }
    return q;
}


static inline short qdigest_alloc(Qdigest *q) {
    short next;
    next = q->free_head;
    q->free_head = q->nodes[q->free_head].next;
    q->num_free--;
    return next;
}

static inline void qdigest_free(Qdigest *q, short index) {
    q->nodes[q->free_tail].next = index;
    q->free_tail = index;
    q->num_free++;
}

//merge two qnode lists, maintaining ascending order
static inline short merge_qnode_lists(Qdigest *q, short a, short b) {
    Qdigest_node *nodes;
    short head, tail, freed;
    nodes = q->nodes;
    if(nodes[a].min < nodes[b].min) { // a is smaller, it should be head
        head = a;
        a = nodes[a].next;
    } else if(nodes[a].min > nodes[b].min) { // b is smaller, it should be head
        head = b;
        b = nodes[b].next;
    } else {
            head = a;
            nodes[a].count += nodes[b].count;
            a = nodes[a].next;
            freed = b;
            b = nodes[b].next;
            qdigest_free(q, freed);
        }
    tail = head;
    while(a && b) {
        if(nodes[a].min < nodes[b].min) {
            nodes[tail].next = a;
            a = nodes[a].next;
        } else if(nodes[a].min > nodes[b].min) {
            nodes[tail].next = b;
            b = nodes[b].next;
        } else {
            nodes[tail].next = a;
            nodes[a].count += nodes[b].count;
            a = nodes[a].next;
            freed = b;
            b = nodes[b].next;
            qdigest_free(q, freed);
        }
        tail = nodes[tail].next;
    }
    //append any leftover nodes to the end of the list
    //if there aren't any leftover nodes, set tail.next to -1
    nodes[tail].next = b == -1 ? a : b;
    return head;
}


static inline void de_duplicate_sorted_list(Qdigest *q, short head) {
    Qdigest_node *nodes;
    short unique_head, scan_head, leftover_head, freed;

    nodes = q->nodes;
    unique_head = head;
    scan_head = nodes[unique_head].next;
    // walk down the list looking for the first pair of identical nodes
    // (this is broken out from the second phase since it is much simpler)
    while(scan_head != -1 && nodes[unique_head].min != nodes[scan_head].min) {
        unique_head = scan_head;
        scan_head = nodes[scan_head].next;
    }
    // now, begin compacting
    while(scan_head != -1) {
        while(nodes[unique_head].min == nodes[scan_head].min) {
            nodes[unique_head].count += nodes[scan_head].count;
            scan_head = nodes[scan_head].next;
        }
        unique_head = nodes[unique_head].next;
        nodes[unique_head].min = nodes[scan_head].min;
        nodes[unique_head].count = nodes[scan_head].count;
    }
    // finally free any left over nodes that were removed by de-duplication
    if(nodes[unique_head].next != -1) {
        leftover_head = nodes[unique_head].next;
        nodes[unique_head].next = -1;
        while(leftover_head != -1) {
            freed = leftover_head;
            leftover_head = nodes[leftover_head].next;
            qdigest_free(q, freed);
        }
    }
}


#define SORTERS_LEN 16

//given an initially unsorted list of qnodes, sort them and merge duplicates
static inline short sort_and_compress(Qdigest *q, short head) {
    Qdigest_node *nodes;
    short unsorted_head, next;
    short sorters[SORTERS_LEN];  
    // sorters[n] is the pointer to a sorted list of length 2 ** n, or -1
    int cur_sorter;
    for(cur_sorter=0; cur_sorter < SORTERS_LEN; cur_sorter++) { sorters[cur_sorter] = -1; }
    nodes = q->nodes;
    unsorted_head = head;
    next = unsorted_head;  // just incase unsorted_head is -1, so next is initialized
    while(unsorted_head != -1) { // combine same length lists until unsorted input exhausted
        if(sorters[0] == -1) {
            sorters[0] = unsorted_head;
            unsorted_head = nodes[unsorted_head].next;
            nodes[sorters[0]].next = -1;
        }
        if(unsorted_head == -1) { break; }
        next = unsorted_head;
        unsorted_head = nodes[unsorted_head].next;
        nodes[next].next = -1;
        for(cur_sorter = 0; cur_sorter < SORTERS_LEN && sorters[cur_sorter] != -1; cur_sorter++) {
            next = merge_qnode_lists(q, next, sorters[cur_sorter]);
        }
        assert(cur_sorter < SORTERS_LEN);
        sorters[cur_sorter] = next;
    }
    // combine differently lengthed lists till there is only one
    // find the first valid linked list in sorters
    for(cur_sorter = 0; cur_sorter < SORTERS_LEN; cur_sorter++) {
        if(sorters[cur_sorter] != -1) {
            next = sorters[cur_sorter];
            break;
        }
    }
    // roll up all the elements in sorters
    for(; cur_sorter < SORTERS_LEN; cur_sorter++) {
        if(sorters[cur_sorter] == -1) { continue; }
        next = merge_qnode_lists(q, sorters[cur_sorter], next);
    }
    // de-duplicate nodes
    de_duplicate_sorted_list(q, next);
    return next;
}


static inline void compress_generation(Qdigest *q, short cur_nodes, short parents, short generation) {
    unsigned int mask;
    short cur, parents_cur, sibling, parent, parents_prev, prev, next;
    unsigned long long count;
    unsigned long long min_count;
    Qdigest_node *nodes;

    nodes = q->nodes;
    min_count = q->n / q->k;
    cur = cur_nodes;
    parents_cur = parents;
    prev = 0;
    parents_prev = 0;
    mask = 0xFFFFFFFF - ((1 << (generation + 2)) - 1);

    while(cur) {
        // test if current node has a sibling
        if(nodes[cur].next && (nodes[cur].min & mask) == (nodes[nodes[cur].next].min & mask)) {
            sibling = nodes[cur].next;
        } else {
            sibling = 0;
        }
        // walk parents to the last one that could possibly be cur parent
        while(parents_cur && (int)(nodes[cur].min & mask) > nodes[parents_cur].min) {
            parents_prev = parents_cur;
            parents_cur = nodes[parents_cur].next;
        }
        // test if current node has a parent
        if(parents_cur && nodes[parents_cur].min == (nodes[cur].min & mask)) {
            parent = parents_cur;
        } else {
            parent = 0;
        }
        next = nodes[sibling ? sibling : cur].next;
        // QDigest count is count(cur) + count(sibling) + count(parent)
        count = nodes[cur].count;
        if(sibling) { count += nodes[sibling].count; }
        if(parent) { count += nodes[parent].count; }
        if(count < min_count) {
            // remove node from current generation list
            if(prev) {
                nodes[prev].next = next;
            } else {
                cur_nodes = next;
            }
            // throw away sibling if one exists
            if(sibling) { qdigest_free(q, sibling); }
            // if there is a parent, throw away cur as well
            if(parent) {
                nodes[parent].count = count;
                qdigest_free(q, cur);
            } else {
                nodes[cur].count = count;
                // reduce granularity to next generation up's
                nodes[cur].min &= mask;
                nodes[cur].next = parents_cur;
                // splice into parents list if possible
                if(parents_prev) {
                    nodes[parents_prev].next = cur;
                } else { // if comes before current parents head, replace
                    parents = cur;
                }
                parents_prev = cur;
            }
            // check if we still need to free more nodes
            if(q->k <= q->num_free) { break; }
        } else { // if nothing merged, advance the prev pointer
            prev = sibling ? sibling : cur;
        }
        cur = next;
    }
    q->generations[generation] = cur_nodes;
    q->generations[generation + 1] = parents;
}

static inline void compress(Qdigest *q) {
    short i;  // don't use char to avoid compiler warnings, even tho max is 32
    short new_nodes;
    new_nodes = sort_and_compress(q, q->new_head);
    q->new_head = q->new_tail = 0;
    q->generations[0] = merge_qnode_lists(
        q, q->generations[0], new_nodes);
    for(i=0; q->num_free < q->k && i < 31; i++) {
        if(q->generations[i] == 0) { continue; }
        compress_generation(q, q->generations[i], q->generations[i+1], i);
    }
}

// this is probably very bad
static union converter {
    unsigned int intval;
    float floatval;
};


void qdigest_update(Qdigest *q, float x) {
    short new_node;
    union converter convert;
    if(!q->num_free) {
        compress(q);
    }
    new_node = qdigest_alloc(q);
    convert.floatval = x;
    q->nodes[new_node].min = convert.intval;
    if(q->new_tail) {
        q->nodes[q->new_tail].next = new_node;
        q->new_tail = new_node;
    } else {
        q->new_head = q->new_tail = new_node;
    }
}


// float must be numeric (not NaN, +inf, or -inf)
// in terms of bits, (X = 1 or 0, doesn't matter)
//  NaN = X-11111111-XXXXXXXXXXXXXXXXXXXXXXX
// +inf = 0-11111111-00000000000000000000000
// -inf = 1-11111111-00000000000000000000000
// this can be reduced to the rule that there are two non-numerical ranges
// 7F80 0000 to 7FFF FFFF, and FF80 0000 to FFFF FFFF
// put another way, 0000 0000 to 7F7F FFFF and 8000 0000 to FF7F FFFF are the valid ranges
static inline float int_bits2numeric_float(unsigned int val) {
    union converter convert;
    if((val & 0x80000000) && val > 0xFF7FFFFF) {
        val = 0xFF7FFFFF;
    } else if(val > 0x7F7FFFFF) {
        val = 0x7F7FFFFF;
    }
    convert.intval = val;
    return convert.floatval;
}

static inline PyObject* _pack_list(Qdigest_node *nodes, short head, int generation) {
    PyObject* ret;
    int size, i;
    short walker;
    if(head == -1) {
        ret = PyTuple_New(0);
    } else {
        // walk once to determine size
        size = 1;
        for(walker = head; walker != -1; walker = nodes[walker].next) {
            size++;
        }
        // walk again to put the values into the tuple
        ret = PyTuple_New(size);
        walker = head;
        for(i=0; i<size; i++) {
            PyTuple_SetItem(ret, (Py_ssize_t)i, 
                PyTuple_Pack(3,
                    PyFloat_FromDouble(int_bits2numeric_float(nodes[walker].min)), 
                    PyFloat_FromDouble(int_bits2numeric_float(nodes[walker].min + (1 << generation) - 1)),
                    PyLong_FromUnsignedLongLong(nodes[walker].count))
            );
            walker = nodes[walker].next;
        }        
    }
    return ret;
}

// NOT a Python function -- a C helper to dump the internal state of the qdigest into python tuples
PyObject* qdigest_dumpstate(Qdigest *q) {
    int i;
    Qdigest_node *nodes;
    PyObject *generations, *in_buff, *ret;
    nodes = q->nodes;
    generations = PyTuple_New(32);
    for(i=0; i < 32; i++) {
        PyTuple_SetItem(generations, (Py_ssize_t)i, _pack_list(nodes, q->generations[i], i));
    }
    in_buff = _pack_list(nodes, q->new_head, 0);
    ret = PyTuple_Pack(2, in_buff, generations);
    return ret;
}
