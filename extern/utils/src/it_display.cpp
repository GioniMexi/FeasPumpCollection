/**
 * @file it_display.cpp
 * @brief Iteration Display Class
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#include "utils/it_display.h"
#include <iostream>
#include <algorithm>

namespace dominiqs {

bool IterationDisplay::addColumn(const std::string& name, int p, int w, int prec, bool v, const std::string& d)
{
	for (auto& c: columns)
	{
		if (c.second->name == name) return false;
	}
	columns[p] = ColumnPtr(new Column(name, p, w, prec, v, d));
	return true;
}

void IterationDisplay::removeColumn(const std::string& name)
{
	auto c = std::find_if(columns.begin(), columns.end(), [&name](const ColumnMap::value_type& itr) { return (itr.second->name == name); });
	columns.erase(c);
}

void IterationDisplay::clear()
{
	columns.clear();
}

void IterationDisplay::setVisible(const std::string& name, bool visible)
{
	for (auto& c: columns)
	{
		if (c.second->name == name)
		{
			c.second->visible = visible;
			break;
		}
	}
}

void IterationDisplay::printHeader(std::ostream& out)
{
	if (columns.empty()) return;
	for (auto& c: columns)
	{
		ColumnPtr col = c.second;
		if (col->visible)
		{
			out << col->formatHeader(col->name);
		}
	}
	out << std::endl;
}

void IterationDisplay::resetIteration()
{
	current.clear();
	marked = false;
}

void IterationDisplay::printIteration(std::ostream& out)
{
	if (columns.empty()) return;
	for (auto& c: columns)
	{
		ColumnPtr col = c.second;
		if (col->visible)
		{
			if (current.find(col->name) != current.end())  out << current[col->name];
			else                                           out << col->formatHeader("");
		}
	}
	out << std::endl;
}

} // namespace dominiqs
